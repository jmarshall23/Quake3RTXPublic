#include "Common.hlsl"

StructuredBuffer<STriVertex> BTriVertex : register(t0);
Texture2D<float4> MegaTexture : register(t1);
StructuredBuffer<SInstanceProperties> BInstanceProperties : register(t2);

float2 CalcUV(int vertId, Attributes attrib)
{
	float3 barycentrics = float3(1.f - attrib.bary.x - attrib.bary.y, attrib.bary.x, attrib.bary.y);

	float u = 0, v = 0;
	for(int i = 0; i < 3; i++)
	{
	u += BTriVertex[vertId + i].st.x * barycentrics[i];
	v += BTriVertex[vertId + i].st.y * barycentrics[i];
	}
	
	u = frac(u) / 16384;
	v = frac(v) / 16384;

	u = u * BTriVertex[vertId + 0].vtinfo.z;
	v = v * BTriVertex[vertId + 0].vtinfo.w;
	
	u = u + (BTriVertex[vertId + 0].vtinfo.x / 16384);
	v = v + (BTriVertex[vertId + 0].vtinfo.y / 16384);
	
	return float2(u, v);
}

[shader("anyhit")] void InteractionAnyHit(inout HitInfo payload, Attributes attrib)
{    
	uint vertId = BInstanceProperties[InstanceID()].startVertex + (3 * PrimitiveIndex());

	float3 worldOrigin = WorldRayOrigin() + RayTCurrent() * WorldRayDirection();

	float u, v;
	float2 result;
	result = CalcUV(vertId, attrib);
	u = result.x;
	v = result.y;
	
	float4 megaColor = MegaTexture.Load(int3(u * 16384, v * 16384, 0)); //normalize(BTriVertex[vertId + 0].vertex) * 4;
	if(megaColor.w != 1)
	{
		payload.colorAndDistance += float4(megaColor.xyz, 1.0);
		IgnoreHit();
		return;
	}
	
	AcceptHitAndEndSearch();
}