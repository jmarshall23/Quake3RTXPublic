// Hit information, aka ray payload
// This sample only carries a shading color and hit distance.
// Note that the payload should be kept as small as possible,
// and that its size must be declared in the corresponding
// D3D12_RAYTRACING_SHADER_CONFIG pipeline subobjet.
struct HitInfo {
  float4 colorAndDistance;
  float4 lightColor;
  float4 worldOrigin;
  float4 worldNormal;
};

// Attributes output by the raytracing when hitting a surface,
// here the barycentric coordinates
struct Attributes {
  float2 bary;
};


struct SecondHitInfo {
   float4 payload_color;
   float4 payload_vert_info; 
};

struct STriVertex {
  float3 vertex;
  float3 st;
  float3 normal;
  float4 vtinfo;
  float4 tangent;
  float4 binormal;
};

struct SInstanceProperties
{
	int startVertex;
};
