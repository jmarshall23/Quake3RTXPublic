// d3d12_raytrace_mesh.cpp
//

#include "d3d12_local.h"
#include "nv_helpers_dx12/BottomLevelASGenerator.h"
#include "nv_helpers_dx12/TopLevelASGenerator.h"
#include <vector>

#define Vector2Subtract(a,b,c)  ((c)[0]=(a)[0]-(b)[0],(c)[1]=(a)[1]-(b)[1])

std::vector<dxrMesh_t *> dxrMeshList;


ComPtr<ID3D12Resource> m_vertexBuffer;
D3D12_VERTEX_BUFFER_VIEW m_vertexBufferView;
std::vector<dxrVertex_t> sceneVertexes;

/*
=============
GL_CalcTangentSpace
Tr3B - recoded from Nvidia's SDK
=============
*/
void GL_CalcTangentSpace(vec3_t tangent, vec3_t binormal, vec3_t normal,
	const vec3_t v0, const vec3_t v1, const vec3_t v2, const vec2_t t0, const vec2_t t1, const vec2_t t2)
{
	vec3_t          cp, e0, e1;
	vec3_t          faceNormal;

	VectorSet(e0, v1[0] - v0[0], t1[0] - t0[0], t1[1] - t0[1]);
	VectorSet(e1, v2[0] - v0[0], t2[0] - t0[0], t2[1] - t0[1]);

	CrossProduct(e0, e1, cp);
	if (fabs(cp[0]) > 10e-6)
	{
		tangent[0] = -cp[1] / cp[0];
		binormal[0] = -cp[2] / cp[0];
	}

	e0[0] = v1[1] - v0[1];
	e1[0] = v2[1] - v0[1];

	CrossProduct(e0, e1, cp);
	if (fabs(cp[0]) > 10e-6)
	{
		tangent[1] = -cp[1] / cp[0];
		binormal[1] = -cp[2] / cp[0];
	}

	e0[0] = v1[2] - v0[2];
	e1[0] = v2[2] - v0[2];

	CrossProduct(e0, e1, cp);
	if (fabs(cp[0]) > 10e-6)
	{
		tangent[2] = -cp[1] / cp[0];
		binormal[2] = -cp[2] / cp[0];
	}

	VectorNormalizeFast(tangent);
	VectorNormalizeFast(binormal);

	// normal...
	// compute the cross product TxB
	CrossProduct(tangent, binormal, normal);
	VectorNormalizeFast(normal);

	// Gram-Schmidt orthogonalization process for B
	// compute the cross product B=NxT to obtain 
	// an orthogonal basis
	CrossProduct(normal, tangent, binormal);

	// compute the face normal based on vertex points
	VectorSubtract(v2, v0, e0);
	VectorSubtract(v1, v0, e1);
	CrossProduct(e0, e1, faceNormal);

	VectorNormalizeFast(faceNormal);

	if (DotProduct(normal, faceNormal) < 0)
	{
		VectorInverse(normal);
		//VectorInverse(tangent);
		//VectorInverse(binormal);
	}
}


/*
=============
GL_RaytraceSurfaceGrid
=============
*/
void GL_RaytraceSurfaceGrid(dxrMesh_t* mesh, msurface_t* fa, srfGridMesh_t* cv) {
	int		i, j;
	//float* xyz;
	//float* texCoords;
	//float* normal;
	unsigned char* color;
	drawVert_t* dv;
	int		rows, irows, vrows;
	int		used;
	int		widthTable[MAX_GRID_SIZE];
	int		heightTable[MAX_GRID_SIZE];
	float	lodError;
	int		lodWidth, lodHeight;
	int		dlightBits;
	int* vDlightBits;
	qboolean	needsNormal;

	dxrSurface_t surf;

	int materialInfo = 0;

	float x, y, w, h;

	if (fa->shader == NULL)
		return;

	x = fa->shader->atlas_x;
	y = fa->shader->atlas_y;
	w = fa->shader->atlas_width;
	h = fa->shader->atlas_height;

	if (!mesh->alphaSurface)
	{
		mesh->alphaSurface = fa->shader->alphaSurface;
	}

	dlightBits = cv->dlightBits[backEnd.smpFrame];
	tess.dlightBits |= dlightBits;

	// determine the allowable discrepance
	lodError = 0;

	// determine which rows and columns of the subdivision
	// we are actually going to use
	widthTable[0] = 0;
	lodWidth = 1;
	for (i = 1; i < cv->width - 1; i++) {
		if (cv->widthLodError[i] <= lodError) {
			widthTable[lodWidth] = i;
			lodWidth++;
		}
	}
	widthTable[lodWidth] = cv->width - 1;
	lodWidth++;

	heightTable[0] = 0;
	lodHeight = 1;
	for (i = 1; i < cv->height - 1; i++) {
		if (cv->heightLodError[i] <= lodError) {
			heightTable[lodHeight] = i;
			lodHeight++;
		}
	}
	heightTable[lodHeight] = cv->height - 1;
	lodHeight++;

	surf.startVertex = mesh->meshVertexes.size();
	surf.numVertexes = 0;

	surf.numIndexes = 0;
	surf.startIndex = mesh->meshIndexes.size();

	// very large grids may have more points or indexes than can be fit
	// in the tess structure, so we may have to issue it in multiple passes

	used = 0;
	rows = 0;
	while (used < lodHeight - 1) {
		// see how many rows of both verts and indexes we can add without overflowing
		do {
			vrows = (SHADER_MAX_VERTEXES - surf.numVertexes) / lodWidth;
			irows = (SHADER_MAX_INDEXES - surf.numIndexes) / (lodWidth * 6);

			// if we don't have enough space for at least one strip, flush the buffer
			if (vrows < 2 || irows < 1) {
				//RB_EndSurface();
				//RB_BeginSurface(tess.shader, tess.fogNum );
			}
			else {
				break;
			}
		} while (1);

		rows = irows;
		if (vrows < irows + 1) {
			rows = vrows - 1;
		}
		if (used + rows > lodHeight) {
			rows = lodHeight - used;
		}

		//xyz = tess.xyz[numVertexes];
		//normal = tess.normal[numVertexes];
		//texCoords = tess.texCoords[numVertexes][0];
		//color = (unsigned char*)&tess.vertexColors[numVertexes];
		//vDlightBits = &tess.vertexDlightBits[numVertexes];
		//needsNormal = tess.shader->needsNormal;

		int startVertex = surf.numVertexes;

		for (i = 0; i < rows; i++) {
			for (j = 0; j < lodWidth; j++) {
				dv = cv->verts + heightTable[used + i] * cv->width
					+ widthTable[j];

				dxrVertex_t v;

				v.xyz[0] = dv->xyz[0];
				v.xyz[1] = dv->xyz[1];
				v.xyz[2] = dv->xyz[2];
				v.st[0] = dv->st[0];
				v.st[1] = dv->st[1];
				v.st[2] = materialInfo;
				v.vtinfo[0] = x;
				v.vtinfo[1] = y;
				v.vtinfo[2] = w;
				v.vtinfo[3] = h;

				mesh->meshVertexes.push_back(v);
				surf.numVertexes++;
			}
		}


		// add the indexes
		{
			int		w, h;

			h = rows - 1;
			w = lodWidth - 1;
			for (i = 0; i < h; i++) {
				for (j = 0; j < w; j++) {
					int		v1, v2, v3, v4;

					// vertex order to be reckognized as tristrips
					v1 = startVertex + i * lodWidth + j + 1;
					v2 = v1 - 1;
					v3 = v2 + lodWidth;
					v4 = v3 + 1;

					mesh->meshIndexes.push_back(surf.startVertex + v2);
					mesh->meshIndexes.push_back(surf.startVertex + v3);
					mesh->meshIndexes.push_back(surf.startVertex + v1);

					mesh->meshIndexes.push_back(surf.startVertex + v1);
					mesh->meshIndexes.push_back(surf.startVertex + v3);
					mesh->meshIndexes.push_back(surf.startVertex + v4);

					surf.numIndexes += 6;
				}
			}
		}

		used += rows - 1;
	}	

	mesh->meshSurfaces.push_back(surf);
}


void GL_LoadBottomLevelAccelStruct(dxrMesh_t* mesh, msurface_t* surfaces, int numSurfaces, int bModelIndex) {
	mesh->startSceneVertex = sceneVertexes.size();

	for (int i = 0; i < numSurfaces; i++)
	{
		msurface_t* fa = &surfaces[i];
		srfTriangles_t* tri = (srfTriangles_t*)fa->data;
		if (tri == NULL) {
			continue;
		}
		if (tri->surfaceType == SF_SKIP) {
			continue;
		}

		if (tri->surfaceType == SF_GRID) {
			GL_RaytraceSurfaceGrid(mesh, fa, (srfGridMesh_t*)fa->data);
			continue;
		}

		dxrSurface_t surf;

		int materialInfo = 1;

		float x, y, w, h;

		if (fa->shader == NULL)
			continue;

		if (fa->shader->defaultShader)
			continue;

		if (strstr(fa->shader->name, "fog")) {
			continue;
		}

		if (strstr(fa->shader->name, "flame")) {
			continue;
		}

		if (strstr(fa->shader->name, "skies") || strstr(fa->shader->name, "light") || strstr(fa->shader->name, "lava")) {
			materialInfo = 2;
		}

		if (strstr(fa->shader->name, "sfx/beam")) {
			continue;
		}

		if (fa->shader->hasRaytracingReflection)
			materialInfo = 5;


		x = fa->shader->atlas_x;
		y = fa->shader->atlas_y;
		w = fa->shader->atlas_width;
		h = fa->shader->atlas_height;

		if (!mesh->alphaSurface)
		{
			mesh->alphaSurface = fa->shader->alphaSurface;
		}

		if (strstr(fa->shader->name, "lavahelldark")) {
			GL_FindMegaTile("lavahell", x, y, w, h);
		}

		//if (w == 0 || h == 0) {
		//	continue;
		//}

		if((fa->bmodelindex > 2 && bModelIndex == -1) || (bModelIndex >= 0 && fa->bmodelindex != bModelIndex)) {
			continue;
		}

		surf.startVertex = mesh->meshVertexes.size();
		surf.numVertexes = 0;
		for (int d = 0; d < tri->numVerts; d++) {
			dxrVertex_t v;

			v.xyz[0] = tri->verts[d].xyz[0];
			v.xyz[1] = tri->verts[d].xyz[1];
			v.xyz[2] = tri->verts[d].xyz[2];
			v.st[0] = tri->verts[d].st[0];
			v.st[1] = tri->verts[d].st[1];
			v.normal[0] = tri->verts[d].normal[0];
			v.normal[1] = tri->verts[d].normal[1];
			v.normal[2] = tri->verts[d].normal[2];
			v.st[2] = materialInfo;
			v.vtinfo[0] = x;
			v.vtinfo[1] = y;
			v.vtinfo[2] = w;
			v.vtinfo[3] = h;

			mesh->meshVertexes.push_back(v);
			surf.numVertexes++;
		}

		surf.numIndexes = 0;
		surf.startIndex = mesh->meshIndexes.size();
		for (int d = 0; d < tri->numIndexes; d++)
		{
			mesh->meshIndexes.push_back(surf.startVertex + tri->indexes[d]);
			surf.numIndexes++;
		}

		mesh->meshSurfaces.push_back(surf);
	}

	// Calculate the normals
	for (int i = 0; i < mesh->meshSurfaces.size(); i++)
	{
		mesh->meshSurfaces[i].startMegaVertex = mesh->meshTriVertexes.size();

		for (int d = 0; d < mesh->meshSurfaces[i].numIndexes; d+=3)
		{			
			int idx0 = mesh->meshIndexes[mesh->meshSurfaces[i].startIndex + d + 0];
			int idx1 = mesh->meshIndexes[mesh->meshSurfaces[i].startIndex + d + 1];
			int idx2 = mesh->meshIndexes[mesh->meshSurfaces[i].startIndex + d + 2];

			dxrVertex_t* vert0 = &mesh->meshVertexes[idx0];
			dxrVertex_t* vert1 = &mesh->meshVertexes[idx1];
			dxrVertex_t* vert2 = &mesh->meshVertexes[idx2];

			vec3_t tangent, binormal, normal;
			GL_CalcTangentSpace(tangent, binormal, normal, vert0->xyz, vert1->xyz, vert2->xyz, vert0->st, vert1->st, vert2->st);

			memcpy(vert0->normal, normal, sizeof(float) * 3);
			memcpy(vert1->normal, normal, sizeof(float) * 3);
			memcpy(vert2->normal, normal, sizeof(float) * 3);

			memcpy(vert0->binormal, binormal, sizeof(float) * 3);
			memcpy(vert1->binormal, binormal, sizeof(float) * 3);
			memcpy(vert2->binormal, binormal, sizeof(float) * 3);

			memcpy(vert0->tangent, tangent, sizeof(float) * 3);
			memcpy(vert1->tangent, tangent, sizeof(float) * 3);
			memcpy(vert2->tangent, tangent, sizeof(float) * 3);
		}
	}

	// TODO: Use a index buffer here : )
	{
		for (int i = 0; i < mesh->meshSurfaces.size(); i++)
		{
			mesh->meshSurfaces[i].startMegaVertex = mesh->meshTriVertexes.size();

			for (int d = 0; d < mesh->meshSurfaces[i].numIndexes; d++)
			{
				int indexId = mesh->meshSurfaces[i].startIndex + d;
				int idx = mesh->meshIndexes[indexId];

				mesh->meshTriVertexes.push_back(mesh->meshVertexes[idx]);
				sceneVertexes.push_back(mesh->meshVertexes[idx]);
				mesh->numSceneVertexes++;
			}
		}
	}	
}

void *GL_LoadDXRMesh(msurface_t *surfaces, int numSurfaces, int bModelIndex)  {
	dxrMesh_t* mesh = new dxrMesh_t();
	
	//mesh->meshId = dxrMeshList.size();
	
	GL_LoadBottomLevelAccelStruct(mesh, surfaces, numSurfaces, bModelIndex);

	dxrMeshList.push_back(mesh);

	return mesh;
}

/*
** LerpMeshVertexes
*/
static void LerpMeshVertexes(md3Surface_t* surf, float backlerp, int frame, int oldframe, dxrVertex_t* vertexes, float x, float y, float w, float h)
{
	short* oldXyz, * newXyz, * oldNormals, * newNormals;
	//float* outXyz, * outNormal;
	float	oldXyzScale, newXyzScale;
	float	oldNormalScale, newNormalScale;
	int		vertNum;
	unsigned lat, lng;
	int		numVerts;
	float* texCoords;

	//outXyz = tess.xyz[tess.numVertexes];
	//outNormal = tess.normal[tess.numVertexes];

	newXyz = (short*)((byte*)surf + surf->ofsXyzNormals)
		+ (frame * surf->numVerts * 4);
	newNormals = newXyz + 3;

	newXyzScale = MD3_XYZ_SCALE * (1.0 - backlerp);
	newNormalScale = 1.0 - backlerp;

	numVerts = surf->numVerts;
	texCoords = (float*)((byte*)surf + surf->ofsSt);
	//
	// just copy the vertexes
	//
	for (vertNum = 0; vertNum < numVerts; vertNum++,
		newXyz += 4, newNormals += 4, vertexes++, texCoords += 2)
	{

		vertexes->xyz[0] = newXyz[0] * newXyzScale;
		vertexes->xyz[1] = newXyz[1] * newXyzScale;
		vertexes->xyz[2] = newXyz[2] * newXyzScale;

		lat = (newNormals[0] >> 8) & 0xff;
		lng = (newNormals[0] & 0xff);
		lat *= (FUNCTABLE_SIZE / 256);
		lng *= (FUNCTABLE_SIZE / 256);

		// decode X as cos( lat ) * sin( long )
		// decode Y as sin( lat ) * sin( long )
		// decode Z as cos( long )

		vertexes->normal[0] = tr.sinTable[(lat + (FUNCTABLE_SIZE / 4)) & FUNCTABLE_MASK] * tr.sinTable[lng];
		vertexes->normal[1] = tr.sinTable[lat] * tr.sinTable[lng];
		vertexes->normal[2] = tr.sinTable[(lng + (FUNCTABLE_SIZE / 4)) & FUNCTABLE_MASK];

		vertexes->st[0] = texCoords[0];
		vertexes->st[1] = texCoords[1];
		vertexes->st[2] = 1;

		vertexes->vtinfo[0] = x;
		vertexes->vtinfo[1] = y;
		vertexes->vtinfo[2] = w;
		vertexes->vtinfo[3] = h;
	}
	/*
		if (backlerp == 0) {
	#if idppc_altivec
			vector signed short newNormalsVec0;
			vector signed short newNormalsVec1;
			vector signed int newNormalsIntVec;
			vector float newNormalsFloatVec;
			vector float newXyzScaleVec;
			vector unsigned char newNormalsLoadPermute;
			vector unsigned char newNormalsStorePermute;
			vector float zero;

			newNormalsStorePermute = vec_lvsl(0, (float*)&newXyzScaleVec);
			newXyzScaleVec = *(vector float*) & newXyzScale;
			newXyzScaleVec = vec_perm(newXyzScaleVec, newXyzScaleVec, newNormalsStorePermute);
			newXyzScaleVec = vec_splat(newXyzScaleVec, 0);
			newNormalsLoadPermute = vec_lvsl(0, newXyz);
			newNormalsStorePermute = vec_lvsr(0, outXyz);
			zero = (vector float)vec_splat_s8(0);
			//
			// just copy the vertexes
			//
			for (vertNum = 0; vertNum < numVerts; vertNum++,
				newXyz += 4, newNormals += 4,
				outXyz += 4, outNormal += 4)
			{
				newNormalsLoadPermute = vec_lvsl(0, newXyz);
				newNormalsStorePermute = vec_lvsr(0, outXyz);
				newNormalsVec0 = vec_ld(0, newXyz);
				newNormalsVec1 = vec_ld(16, newXyz);
				newNormalsVec0 = vec_perm(newNormalsVec0, newNormalsVec1, newNormalsLoadPermute);
				newNormalsIntVec = vec_unpackh(newNormalsVec0);
				newNormalsFloatVec = vec_ctf(newNormalsIntVec, 0);
				newNormalsFloatVec = vec_madd(newNormalsFloatVec, newXyzScaleVec, zero);
				newNormalsFloatVec = vec_perm(newNormalsFloatVec, newNormalsFloatVec, newNormalsStorePermute);
				//outXyz[0] = newXyz[0] * newXyzScale;
				//outXyz[1] = newXyz[1] * newXyzScale;
				//outXyz[2] = newXyz[2] * newXyzScale;

				lat = (newNormals[0] >> 8) & 0xff;
				lng = (newNormals[0] & 0xff);
				lat *= (FUNCTABLE_SIZE / 256);
				lng *= (FUNCTABLE_SIZE / 256);

				// decode X as cos( lat ) * sin( long )
				// decode Y as sin( lat ) * sin( long )
				// decode Z as cos( long )

				outNormal[0] = tr.sinTable[(lat + (FUNCTABLE_SIZE / 4)) & FUNCTABLE_MASK] * tr.sinTable[lng];
				outNormal[1] = tr.sinTable[lat] * tr.sinTable[lng];
				outNormal[2] = tr.sinTable[(lng + (FUNCTABLE_SIZE / 4)) & FUNCTABLE_MASK];

				vec_ste(newNormalsFloatVec, 0, outXyz);
				vec_ste(newNormalsFloatVec, 4, outXyz);
				vec_ste(newNormalsFloatVec, 8, outXyz);
			}

	#else
			//
			// just copy the vertexes
			//
			for (vertNum = 0; vertNum < numVerts; vertNum++,
				newXyz += 4, newNormals += 4, vertexes++)
			{

				vertexes->xyz[0] = newXyz[0] * newXyzScale;
				vertexes->xyz[1] = newXyz[1] * newXyzScale;
				vertexes->xyz[2] = newXyz[2] * newXyzScale;

				lat = (newNormals[0] >> 8) & 0xff;
				lng = (newNormals[0] & 0xff);
				lat *= (FUNCTABLE_SIZE / 256);
				lng *= (FUNCTABLE_SIZE / 256);

				// decode X as cos( lat ) * sin( long )
				// decode Y as sin( lat ) * sin( long )
				// decode Z as cos( long )

				outNormal[0] = tr.sinTable[(lat + (FUNCTABLE_SIZE / 4)) & FUNCTABLE_MASK] * tr.sinTable[lng];
				outNormal[1] = tr.sinTable[lat] * tr.sinTable[lng];
				outNormal[2] = tr.sinTable[(lng + (FUNCTABLE_SIZE / 4)) & FUNCTABLE_MASK];
			}
	#endif
		}
		else {
			//
			// interpolate and copy the vertex and normal
			//
			oldXyz = (short*)((byte*)surf + surf->ofsXyzNormals)
				+ (oldframe * surf->numVerts * 4);
			oldNormals = oldXyz + 3;

			oldXyzScale = MD3_XYZ_SCALE * backlerp;
			oldNormalScale = backlerp;

			for (vertNum = 0; vertNum < numVerts; vertNum++,
				oldXyz += 4, newXyz += 4, oldNormals += 4, newNormals += 4,
				outXyz += 4, outNormal += 4)
			{
				vec3_t uncompressedOldNormal, uncompressedNewNormal;

				// interpolate the xyz
				outXyz[0] = oldXyz[0] * oldXyzScale + newXyz[0] * newXyzScale;
				outXyz[1] = oldXyz[1] * oldXyzScale + newXyz[1] * newXyzScale;
				outXyz[2] = oldXyz[2] * oldXyzScale + newXyz[2] * newXyzScale;

				// FIXME: interpolate lat/long instead?
				lat = (newNormals[0] >> 8) & 0xff;
				lng = (newNormals[0] & 0xff);
				lat *= 4;
				lng *= 4;
				uncompressedNewNormal[0] = tr.sinTable[(lat + (FUNCTABLE_SIZE / 4)) & FUNCTABLE_MASK] * tr.sinTable[lng];
				uncompressedNewNormal[1] = tr.sinTable[lat] * tr.sinTable[lng];
				uncompressedNewNormal[2] = tr.sinTable[(lng + (FUNCTABLE_SIZE / 4)) & FUNCTABLE_MASK];

				lat = (oldNormals[0] >> 8) & 0xff;
				lng = (oldNormals[0] & 0xff);
				lat *= 4;
				lng *= 4;

				uncompressedOldNormal[0] = tr.sinTable[(lat + (FUNCTABLE_SIZE / 4)) & FUNCTABLE_MASK] * tr.sinTable[lng];
				uncompressedOldNormal[1] = tr.sinTable[lat] * tr.sinTable[lng];
				uncompressedOldNormal[2] = tr.sinTable[(lng + (FUNCTABLE_SIZE / 4)) & FUNCTABLE_MASK];

				outNormal[0] = uncompressedOldNormal[0] * oldNormalScale + uncompressedNewNormal[0] * newNormalScale;
				outNormal[1] = uncompressedOldNormal[1] * oldNormalScale + uncompressedNewNormal[1] * newNormalScale;
				outNormal[2] = uncompressedOldNormal[2] * oldNormalScale + uncompressedNewNormal[2] * newNormalScale;

				//			VectorNormalize (outNormal);
			}
			VectorArrayNormalize((vec4_t*)tess.normal[tess.numVertexes], numVerts);
		}
	*/
}

void* GL_LoadMD3RaytracedMesh(md3Header_t* mod, int frame) {
	dxrMesh_t* mesh = new dxrMesh_t();

	mesh->alphaSurface = qtrue;

	//mesh->meshId = dxrMeshList.size();
	mesh->startSceneVertex = sceneVertexes.size();
	mesh->numSceneVertexes = 0;

	float x, y, w, h;
	char textureName[512];

	md3Surface_t* surf = (md3Surface_t*)((byte*)mod + mod->ofsSurfaces);

	for (int i = 0; i < mod->numSurfaces; i++) {
		md3Shader_t* shader = (md3Shader_t*)((byte*)surf + surf->ofsShaders);

		//if (!mesh->alphaSurface)
		//{
		//	mesh->alphaSurface = shader->alphaSurface;
		//}

		COM_StripExtension(COM_SkipPath((char*)shader->name), textureName);
		GL_FindMegaTile(textureName, &x, &y, &w, &h);

		int startVert = mesh->meshVertexes.size();
		dxrVertex_t* meshVertexes = new dxrVertex_t[surf->numVerts];

		LerpMeshVertexes(surf, 0.0f, frame, frame, meshVertexes, x, y, w, h);

		int indexes = surf->numTriangles * 3;
		int* triangles = (int*)((byte*)surf + surf->ofsTriangles);

		for (int j = 0; j < indexes; j++) {
			int tri = triangles[j];
			dxrVertex_t v = meshVertexes[tri];

			mesh->meshTriVertexes.push_back(v);
			sceneVertexes.push_back(v);
			mesh->numSceneVertexes++;
		}

		delete[] meshVertexes;

		// find the next surface
		surf = (md3Surface_t*)((byte*)surf + surf->ofsEnd);
	}

	// Calculate the normals
	{
		for (int i = 0; i < mesh->numSceneVertexes; i += 3)
		{
			float* pA = &sceneVertexes[mesh->startSceneVertex + i + 0].xyz[0];
			float* pC = &sceneVertexes[mesh->startSceneVertex + i + 1].xyz[0];
			float* pB = &sceneVertexes[mesh->startSceneVertex + i + 2].xyz[0];

			float* tA = &sceneVertexes[mesh->startSceneVertex + i + 0].st[0];
			float* tC = &sceneVertexes[mesh->startSceneVertex + i + 1].st[0];
			float* tB = &sceneVertexes[mesh->startSceneVertex + i + 2].st[0];

			vec3_t tangent, binormal, normal;
			GL_CalcTangentSpace(tangent, binormal, normal, sceneVertexes[mesh->startSceneVertex + i + 0].xyz, sceneVertexes[mesh->startSceneVertex + i + 1].xyz, sceneVertexes[mesh->startSceneVertex + i + 2].xyz,
				sceneVertexes[mesh->startSceneVertex + i + 0].st, sceneVertexes[mesh->startSceneVertex + i + 1].st, sceneVertexes[mesh->startSceneVertex + i + 2].st);

			memcpy(sceneVertexes[mesh->startSceneVertex + i + 0].normal, normal, sizeof(float) * 3);
			memcpy(sceneVertexes[mesh->startSceneVertex + i + 1].normal, normal, sizeof(float) * 3);
			memcpy(sceneVertexes[mesh->startSceneVertex + i + 2].normal, normal, sizeof(float) * 3);

			memcpy(sceneVertexes[mesh->startSceneVertex + i + 0].binormal, binormal, sizeof(float) * 3);
			memcpy(sceneVertexes[mesh->startSceneVertex + i + 1].binormal, binormal, sizeof(float) * 3);
			memcpy(sceneVertexes[mesh->startSceneVertex + i + 2].binormal, binormal, sizeof(float) * 3);

			memcpy(sceneVertexes[mesh->startSceneVertex + i + 0].tangent, tangent, sizeof(float) * 3);
			memcpy(sceneVertexes[mesh->startSceneVertex + i + 1].tangent, tangent, sizeof(float) * 3);
			memcpy(sceneVertexes[mesh->startSceneVertex + i + 2].tangent, tangent, sizeof(float) * 3);
		}
	}

	dxrMeshList.push_back(mesh);

	return mesh;
}

void* GL_LoadPolyRaytracedMesh(shader_t *shader, polyVert_t* verts, int numVertexes) {
	dxrMesh_t* mesh = new dxrMesh_t();

	mesh->alphaSurface = qtrue;

	//mesh->meshId = dxrMeshList.size();
	mesh->startSceneVertex = sceneVertexes.size();
	mesh->numSceneVertexes = 0;

	float x, y, w, h;
	char textureName[512];
	COM_StripExtension(COM_SkipPath((char*)shader->name), textureName);
	GL_FindMegaTile(textureName, &x, &y, &w, &h);


	for (int j = 0; j < numVertexes; j++) {
		dxrVertex_t v;

		v.xyz[0] = verts[j].xyz[0];
		v.xyz[1] = verts[j].xyz[1];
		v.xyz[2] = verts[j].xyz[2];
		v.st[0] = verts[j].st[0];
		v.st[1] = verts[j].st[1];
		v.st[2] = 1;
		v.vtinfo[0] = x;
		v.vtinfo[1] = y;
		v.vtinfo[2] = w;
		v.vtinfo[3] = h;

		mesh->meshTriVertexes.push_back(v);
		sceneVertexes.push_back(v);
		mesh->numSceneVertexes++;
	}

	// Calculate the normals
	{
		for (int i = 0; i < mesh->numSceneVertexes; i += 3)
		{
			float* pA = &sceneVertexes[mesh->startSceneVertex + i + 0].xyz[0];
			float* pC = &sceneVertexes[mesh->startSceneVertex + i + 1].xyz[0];
			float* pB = &sceneVertexes[mesh->startSceneVertex + i + 2].xyz[0];

			float* tA = &sceneVertexes[mesh->startSceneVertex + i + 0].st[0];
			float* tC = &sceneVertexes[mesh->startSceneVertex + i + 1].st[0];
			float* tB = &sceneVertexes[mesh->startSceneVertex + i + 2].st[0];

			vec3_t tangent, binormal, normal;
			GL_CalcTangentSpace(tangent, binormal, normal, sceneVertexes[mesh->startSceneVertex + i + 0].xyz, sceneVertexes[mesh->startSceneVertex + i + 1].xyz, sceneVertexes[mesh->startSceneVertex + i + 2].xyz,
				sceneVertexes[mesh->startSceneVertex + i + 0].st, sceneVertexes[mesh->startSceneVertex + i + 1].st, sceneVertexes[mesh->startSceneVertex + i + 2].st);

			memcpy(sceneVertexes[mesh->startSceneVertex + i + 0].normal, normal, sizeof(float) * 3);
			memcpy(sceneVertexes[mesh->startSceneVertex + i + 1].normal, normal, sizeof(float) * 3);
			memcpy(sceneVertexes[mesh->startSceneVertex + i + 2].normal, normal, sizeof(float) * 3);

			memcpy(sceneVertexes[mesh->startSceneVertex + i + 0].binormal, binormal, sizeof(float) * 3);
			memcpy(sceneVertexes[mesh->startSceneVertex + i + 1].binormal, binormal, sizeof(float) * 3);
			memcpy(sceneVertexes[mesh->startSceneVertex + i + 2].binormal, binormal, sizeof(float) * 3);

			memcpy(sceneVertexes[mesh->startSceneVertex + i + 0].tangent, tangent, sizeof(float) * 3);
			memcpy(sceneVertexes[mesh->startSceneVertex + i + 1].tangent, tangent, sizeof(float) * 3);
			memcpy(sceneVertexes[mesh->startSceneVertex + i + 2].tangent, tangent, sizeof(float) * 3);
		}
	}

	dxrMeshList.push_back(mesh);

	return mesh;
}

void GL_FinishVertexBufferAllocation(void) {
//	ThrowIfFailed(m_commandList->Reset(m_commandAllocator.Get(), m_pipelineState.Get()));

	// Create the vertex buffer.
	{
		const UINT vertexBufferSize = sizeof(dxrVertex_t) * sceneVertexes.size();

		// Note: using upload heaps to transfer static data like vert buffers is not 
		// recommended. Every time the GPU needs it, the upload heap will be marshalled 
		// over. Please read up on Default Heap usage. An upload heap is used here for 
		// code simplicity and because there are very few verts to actually transfer.
		ThrowIfFailed(m_device->CreateCommittedResource(
			&CD3DX12_HEAP_PROPERTIES(D3D12_HEAP_TYPE_UPLOAD),
			D3D12_HEAP_FLAG_NONE,
			&CD3DX12_RESOURCE_DESC::Buffer(vertexBufferSize),
			D3D12_RESOURCE_STATE_GENERIC_READ,
			nullptr,
			IID_PPV_ARGS(&m_vertexBuffer)));

		// Copy the triangle data to the vertex buffer.
		UINT8* pVertexDataBegin;
		CD3DX12_RANGE readRange(0, 0);		// We do not intend to read from this resource on the CPU.
		ThrowIfFailed(m_vertexBuffer->Map(0, &readRange, reinterpret_cast<void**>(&pVertexDataBegin)));
		memcpy(pVertexDataBegin, &sceneVertexes[0], sizeof(dxrVertex_t) * sceneVertexes.size());
		m_vertexBuffer->Unmap(0, nullptr);

		// Initialize the vertex buffer view.
		m_vertexBufferView.BufferLocation = m_vertexBuffer->GetGPUVirtualAddress();
		m_vertexBufferView.StrideInBytes = sizeof(dxrVertex_t);
		m_vertexBufferView.SizeInBytes = vertexBufferSize;
	}

	for(int i = 0; i < dxrMeshList.size(); i++)
	{
		dxrMesh_t* mesh = dxrMeshList[i];

		nv_helpers_dx12::BottomLevelASGenerator bottomLevelAS;
		bottomLevelAS.AddVertexBuffer(m_vertexBuffer.Get(), mesh->startSceneVertex * sizeof(dxrVertex_t), mesh->numSceneVertexes, sizeof(dxrVertex_t), NULL, 0, !mesh->alphaSurface);

		// Adding all vertex buffers and not transforming their position.
		//for (const auto& buffer : vVertexBuffers) {
		//	bottomLevelAS.AddVertexBuffer(buffer.first.Get(), 0, buffer.second,
		//		sizeof(Vertex), 0, 0);
		//}

		// The AS build requires some scratch space to store temporary information.
		// The amount of scratch memory is dependent on the scene complexity.
		UINT64 scratchSizeInBytes = 0;
		// The final AS also needs to be stored in addition to the existing vertex
		// buffers. It size is also dependent on the scene complexity.
		UINT64 resultSizeInBytes = 0;

		bottomLevelAS.ComputeASBufferSizes(m_device.Get(), false, &scratchSizeInBytes,
			&resultSizeInBytes);

		// Once the sizes are obtained, the application is responsible for allocating
		// the necessary buffers. Since the entire generation will be done on the GPU,
		// we can directly allocate those on the default heap	
		mesh->buffers.pScratch = nv_helpers_dx12::CreateBuffer(m_device.Get(), scratchSizeInBytes, D3D12_RESOURCE_FLAG_ALLOW_UNORDERED_ACCESS, D3D12_RESOURCE_STATE_COMMON, nv_helpers_dx12::kDefaultHeapProps);
		mesh->buffers.pResult = nv_helpers_dx12::CreateBuffer(m_device.Get(), resultSizeInBytes, D3D12_RESOURCE_FLAG_ALLOW_UNORDERED_ACCESS, D3D12_RESOURCE_STATE_RAYTRACING_ACCELERATION_STRUCTURE, nv_helpers_dx12::kDefaultHeapProps);

		// Build the acceleration structure. Note that this call integrates a barrier
		// on the generated AS, so that it can be used to compute a top-level AS right
		// after this method.

		bottomLevelAS.Generate(m_commandList.Get(), mesh->buffers.pScratch.Get(), mesh->buffers.pResult.Get(), false, nullptr);
	}

	// Flush the command list and wait for it to finish
	//m_commandList->Close();
	//ID3D12CommandList* ppCommandLists[] = { m_commandList.Get() };
	//m_commandQueue->ExecuteCommandLists(1, ppCommandLists);
	//m_fenceValue++;
	//m_commandQueue->Signal(m_fence.Get(), m_fenceValue);
	//
	//m_fence->SetEventOnCompletion(m_fenceValue, m_fenceEvent);
	//WaitForSingleObject(m_fenceEvent, INFINITE);
}