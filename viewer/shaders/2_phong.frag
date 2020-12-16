#version 410

uniform float lightIntensity;
uniform bool blinnPhong;
uniform float shininess;
uniform float eta;
uniform sampler2D shadowMap;

in vec4 eyeVector;
in vec4 lightVector;
in vec4 vertColor;
in vec4 vertNormal;
in vec4 lightSpace;

out vec4 fragColor;

/* ----- Fresnel coefficient functions ----- 
There are many ways to compute the fresnel coefficient:
     - It can be approximated
     - It can be computed for dielectric materials with eta real
     - It can be computed for dielectric and conductors with eta complex
     - eta can be expressed as eta_out / eta_in or eta_in and eta_out can be directly given as inputs of the function
     - eta_in = 1 with air/other interface lead to simplifications
     - It can be a simple float or a vec3 depending if we want to apply the same coefficient to each color channel or not

In the following sections, few functions are given for computing fresnel coefficients, the code and formulas are inspired from : 
     - The class / TP notes
     - This website : https://seblagarde.wordpress.com/2013/04/29/memo-on-fresnel-equations/#more-1921
*/

/*
This function is inspired from the class notes. 
It return a single float and take a single real eta as a parameter.
Here, for the function to work, eta must be eta_out/eta_in
*/
float fresnel_real(vec4 half, vec4 outVect, float eta) {
	float cosTheta = dot(half, outVect);
	float sinTheta = length(cross(half.xyz, outVect.xyz));

	float ci = sqrt(max(eta*eta - sinTheta*sinTheta, 0.0));
	float Fs = pow((cosTheta - ci) / (cosTheta + ci), 2);
	float Fp = pow((eta*eta*cosTheta - ci) / (eta*eta*cosTheta + ci), 2);
	
     float F = (Fs + Fp) / 2.0;
	return F;
}
// -----------------------------------------

vec4 getAmbiant(float k, vec4 C, float I) {
     return k * C * I;
}

vec4 getDiffuse(float k, vec4 C, float I, vec4 n, vec4 L) {
     return k * C * max(dot(n, L), 0.0) * I;
}

/*
Return the specular component for Blinn-Phong shading with a fresnel coefficient with a real eta.
*/
vec4 getSpecular_EtaReal_Phong_Fresnel(vec4 C, float I, vec4 n, vec4 L, vec4 V, float eta_in, float eta_out, float exponent) {
     float eta = eta_out / eta_in;
     vec4 H = normalize(V + L);
     float F = fresnel_real(H, L, eta);

     return F * C * pow(max(dot(n, H), 0.0), exponent) * I;
}

void main( void )
{
     vec4 C = vertColor;

     vec4 n = normalize(vertNormal);
     vec4 L = normalize(lightVector);
     vec4 V = normalize(eyeVector);

     // Ambiant 
     float kA = 0.2;
     vec4 ambiant = getAmbiant(kA, C, lightIntensity);
     
     // Diffuse
     float kD = 0.8;
     vec4 diffuse = getDiffuse(kD, C, lightIntensity, n, L);
     
     // Specular
     vec4 specular = getSpecular_EtaReal_Phong_Fresnel(C, lightIntensity, n, L, V, 1.0, eta, shininess);
     
     vec4 finalColor = ambiant + diffuse + specular;
     fragColor = finalColor;
}
