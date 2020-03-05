#version 330

in vec3 vViewSpacePosition;
in vec3 vViewSpaceNormal;
in vec2 vTexCoords;

uniform vec3 uLightDirection;
uniform vec3 uLightIntensity;

uniform vec4 uBaseColorFactor;


uniform sampler2D uMetallicRoughnessTexture;
uniform float uMetallicFactor;
uniform float uRoughnessFactor;

uniform sampler2D uBaseColorTexture;

uniform sampler2D uEmissiveTexture;
uniform vec3 uEmissiveFactor;

out vec3 fColor;

// Constants
const float GAMMA = 2.2;
const float INV_GAMMA = 1. / GAMMA;
const float M_PI = 3.141592653589793;
const float M_1_PI = 1.0 / M_PI;

// We need some simple tone mapping functions
// Basic gamma = 2.2 implementation
// Stolen here: https://github.com/KhronosGroup/glTF-Sample-Viewer/blob/master/src/shaders/tonemapping.glsl

// linear to sRGB approximation
// see http://chilliant.blogspot.com/2012/08/srgb-approximations-for-hlsl.html
vec3 LINEARtoSRGB(vec3 color)
{
  return pow(color, vec3(INV_GAMMA));
}

// sRGB to linear approximation
// see http://chilliant.blogspot.com/2012/08/srgb-approximations-for-hlsl.html
vec4 SRGBtoLINEAR(vec4 srgbIn)
{
  return vec4(pow(srgbIn.xyz, vec3(GAMMA)), srgbIn.w);
}

void main()
{
	vec4 metallicRoughnessTexture = texture(uMetallicRoughnessTexture, vTexCoords);
	float metallicFactor = uMetallicFactor * metallicRoughnessTexture.b;
	float roughnessFactor = uRoughnessFactor *  metallicRoughnessTexture.g;

	vec3 N = normalize(vViewSpaceNormal);
	vec3 L = uLightDirection;

	vec3 V = normalize(-vViewSpacePosition);
	vec3 H = normalize(L + V);
	float a = roughnessFactor * roughnessFactor;

	vec4 baseColorFromTexture = SRGBtoLINEAR(texture(uBaseColorTexture, vTexCoords));
	vec4 baseColor = uBaseColorFactor * baseColorFromTexture;

	const vec3 dielectricSpecular = vec3(0.04, 0.04, 0.04);
	const vec3 black = vec3(0, 0, 0);
	vec3 F0 = mix(dielectricSpecular, baseColor.rgb, metallicFactor);

	float NdotL = clamp(dot(N, L), 0, 1);
	float NdotV = clamp(dot(N, V), 0, 1);
	float NdotH = clamp(dot(N, H), 0, 1);
	float VdotH = clamp(dot(V, H), 0, 1);

	float baseShlickFactor = 1 - VdotH; //need to be exponentiate by 5
	float shlickFactor = baseShlickFactor * baseShlickFactor; //pow2
	shlickFactor *= shlickFactor; //pow4
	shlickFactor *= baseShlickFactor; //pow5


	vec3 F = F0 + (1 - F0) * shlickFactor;

	float sqrtV = sqrt(NdotV * NdotV * (1 - a*a) + a*a);
	float sqrtL = sqrt(NdotL * NdotL * (1 - a*a) + a*a);
	float Vis = 0.5 / (NdotL * sqrtV + NdotV * sqrtL);

	float D = a*a / (M_1_PI * ( pow( NdotH * NdotH * (a*a-1) +1, 2)));

	vec3 f_diffuse = baseColorFromTexture.rgb * M_1_PI * NdotL;
	vec3 f_specular = F * Vis * D;


	//emissive component
	vec3 emissiveTexture = SRGBtoLINEAR(texture(uEmissiveTexture, vTexCoords)).rgb;
	vec3 emissiveComponent = emissiveTexture * uEmissiveFactor;

	fColor = LINEARtoSRGB((f_diffuse + f_specular) * uLightIntensity * NdotL + emissiveComponent);//vec3(vTexCoords, 0); //for debug
}