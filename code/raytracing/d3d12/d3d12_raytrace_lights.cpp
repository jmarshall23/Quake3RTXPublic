// d3d12_raytrace_lights.cpp
//

#include "d3d12_local.h"
#include "../math/vectormath.h"

#define MAX_DRAW_LIGHTS MAX_DLIGHTS
#define MAX_WORLD_LIGHTS 512

struct glLight_t {
	vec4_t origin_radius;
	vec4_t light_color;
	vec4_t light_normal;
	vec3_t absmin;
	vec3_t absmax;
	refEntity_t* ent;
	int distance;
};

struct sceneLightInfo_t {
	vec4_t origin_radius;
	vec4_t light_color;
	vec4_t light_normal;
};

glLight_t worldLights[MAX_WORLD_LIGHTS];
int numWorldLights = 0;

sceneLightInfo_t *sceneLights = NULL;
tr_buffer* sceneLightInfoBuffer;

/*
===============
GL_ClearLights
===============
*/
void GL_ClearLights(void) {
	memset(&worldLights[0], 0, sizeof(worldLights));
	numWorldLights = 0;
}


void GL_RegisterWorldLight(refEntity_t* ent, float x, float y, float z, float radius, float r, float g, float b, vec3_t lightNormal) {
	glLight_t light;
	light.origin_radius[0] = x;
	light.origin_radius[1] = y;
	light.origin_radius[2] = z;
	light.origin_radius[3] = radius;
	light.absmin[0] = x;
	light.absmin[1] = y;
	light.absmin[2] = z;

	light.absmax[0] = x;
	light.absmax[1] = y;
	light.absmax[2] = z;

	light.ent = ent;

	light.light_color[0] = r;
	light.light_color[1] = g;
	light.light_color[2] = b;

	light.light_normal[0] = lightNormal[0];
	light.light_normal[1] = lightNormal[1];
	light.light_normal[2] = lightNormal[2];

	worldLights[numWorldLights++] = light;
}

void GL_InitLightInfoBuffer(D3D12_CPU_DESCRIPTOR_HANDLE& srvPtr) {
	tr_create_uniform_buffer(renderer, sizeof(sceneLightInfo_t) * MAX_DRAW_LIGHTS, true, &sceneLightInfoBuffer);
	sceneLights = (sceneLightInfo_t *)sceneLightInfoBuffer->cpu_mapped_address;

	D3D12_SHADER_RESOURCE_VIEW_DESC srvDesc = {};
	
	uint32_t bufferSize = ROUND_UP(sizeof(sceneLightInfo_t), D3D12_CONSTANT_BUFFER_DATA_PLACEMENT_ALIGNMENT);
	
	srvDesc.Shader4ComponentMapping = D3D12_DEFAULT_SHADER_4_COMPONENT_MAPPING;
	srvDesc.Format = DXGI_FORMAT_UNKNOWN;
	srvDesc.ViewDimension = D3D12_SRV_DIMENSION_BUFFER;
	srvDesc.Buffer.FirstElement = 0;
	srvDesc.Buffer.NumElements = MAX_DRAW_LIGHTS;
	srvDesc.Buffer.StructureByteStride = sizeof(sceneLightInfo_t);
	srvDesc.Buffer.Flags = D3D12_BUFFER_SRV_FLAG_NONE;
	// Write the per-instance properties buffer view in the heap
	m_device->CreateShaderResourceView(sceneLightInfoBuffer->dx_resource, &srvDesc, srvPtr);
}

int lightSort(const void* a, const void* b) {
	glLight_t light1 = *((glLight_t*)a);
	glLight_t light2 = *((glLight_t*)b);

	return light1.distance - light2.distance;
}

void GL_BuildLightList(float x, float y, float z) {
	int numVisLights = 0;

	memset(sceneLights, 0, sizeof(sceneLightInfo_t) * MAX_DRAW_LIGHTS);
	
	for (int i = 0; i < numWorldLights; i++) {
		glLight_t* ent = &worldLights[i];
		vec3_t viewpos = { x, y, z };

		ent->distance = Distance(ent->origin_radius, viewpos);
	}

	qsort(worldLights, numWorldLights, sizeof(glLight_t), lightSort);

	for(int i = 0; i < numWorldLights; i++) {
		if(numVisLights >= MAX_DRAW_LIGHTS) {
			//Com_Printf("MAX_DRAW_LIGHTS!\n");
			break;
		}
	
		//glLight_t * ent = &worldLights[i];
		//vec3_t viewpos = { x, y, z };
		//byte* pvs = SV_FatPVS(viewpos, cl.worldmodel);
		//
		//int d;
		//for (d = 0; d < ent->num_leafs; d++)
		//	if (pvs[ent->leafnums[d] >> 3] & (1 << (ent->leafnums[d] & 7)))
		//		break;
		//
		//if (d == ent->num_leafs)
		//	continue;
	
		sceneLights[numVisLights].origin_radius[0] = worldLights[i].origin_radius[0];
		sceneLights[numVisLights].origin_radius[1] = worldLights[i].origin_radius[1];
		sceneLights[numVisLights].origin_radius[2] = worldLights[i].origin_radius[2];
		sceneLights[numVisLights].origin_radius[3] = worldLights[i].origin_radius[3];
	
		sceneLights[numVisLights].light_color[0] = worldLights[i].light_color[0];
		sceneLights[numVisLights].light_color[1] = worldLights[i].light_color[1];
		sceneLights[numVisLights].light_color[2] = worldLights[i].light_color[2];

		sceneLights[numVisLights].light_normal[0] = worldLights[i].light_normal[0];
		sceneLights[numVisLights].light_normal[1] = worldLights[i].light_normal[1];
		sceneLights[numVisLights].light_normal[2] = worldLights[i].light_normal[2];

		numVisLights++;
	}
}