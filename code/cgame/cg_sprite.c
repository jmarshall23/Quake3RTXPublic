// cg_sprite.c
//

#include "cg_local.h"

polyVert_t spriteVerts[] = {
	{   0.0,-1.0, -1.0,  0.0, 0.0, 1, 1, 1 , 1},
	{   0.0,1.0,  -1.0,  1.0, 0.0, 1, 1, 1 , 1},
	{   0.0,-1.0,  1.0,  0.0, 1.0, 1, 1, 1 , 1},

	{   0.0,-1.0, 1.0,   0.0, 1.0, 1, 1, 1 , 1},
	{   0.0,1.0, -1.0,   1.0, 0.0, 1, 1, 1 , 1},
	{   0.0,1.0,  1.0,   1.0, 1.0, 1, 1, 1 , 1},
};

/*
===============
CG_InitSmokePuffSprite
===============
*/
void CG_InitSmokePuffSprite(void) {
	qhandle_t shader;
	shader = trap_R_RegisterShader("smokePuff");
	cgs.media.smokePuffModel = trap_R_RegisterCustomModel("_smokePuffMesh", shader, &spriteVerts, 6);
}

/*
===============
CG_InitPlasmaSprite
===============
*/
void CG_InitPlasmaSprite(void) {
	qhandle_t shader;
	shader = trap_R_RegisterShader("sprites/plasma1");	
	cgs.media.plasmaBallModel =  trap_R_RegisterCustomModel("_plasmaSpriteMesh", shader, &spriteVerts, 6);
}

/*
===============
CG_InitSprites
===============
*/
void CG_InitSprites(void) {
	CG_InitPlasmaSprite();
	CG_InitSmokePuffSprite();
}

/*
=================
CG_SpawnSprite
=================
*/
void CG_SpawnSprite(vec3_t origin, float scale, qhandle_t model) {
	refEntity_t	 ent;

	memset(&ent, 0, sizeof(refEntity_t));
	VectorCopy(origin, ent.origin);
	VectorCopy(origin, ent.oldorigin);

	ent.reType = RT_MODEL;
	ent.hModel = model;
	ent.scale = scale;	

	AnglesToAxis(cg.refdefViewAngles, ent.axis);

	trap_R_AddRefEntityToScene(&ent);
}

/*
=================
CG_SpawnSprite
=================
*/
void CG_SpawnSpriteEx(refEntity_t *ent, vec3_t origin, float scale, qhandle_t model) {
	VectorCopy(origin, ent->origin);
	VectorCopy(origin, ent->oldorigin);

	ent->reType = RT_MODEL;
	ent->hModel = model;
	ent->scale = scale;
	AnglesToAxis(cg.refdefViewAngles, ent->axis);	
}