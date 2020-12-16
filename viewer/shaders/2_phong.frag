#version 410

uniform float lightIntensity;
uniform bool blinnPhong;
uniform float shininess;
uniform float eta;
uniform float etak; // Imaginary part when eta is complex
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

In the following sections, few functions are given for computing fresnel coefficients, the code and formulas are inspired from: 
     - The class / TP notes
     - This website : https://seblagarde.wordpress.com/2013/04/29/memo-on-fresnel-equations/#more-1921

Note: The eta_in, eta_out notation might be confusing, eta_out is the index of refraction of the material we are entering/bouncing of off and eta_in is the index of refraction of the medium.
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

/*
This function is taken from https://seblagarde.wordpress.com/2013/04/29/memo-on-fresnel-equations/#more-1921
It return a single float and take a single complex eta (separated in eta the real part and etak the imaginary part) and cosTheta as parameters.
Here, for the function to work, eta + i*etak must be eta_out/eta_in + i*(etak_out/eta_in) 
*/
float fresnel_dieletric_conductor(float eta, float etak, float cosTheta)
{  
   float CosTheta2 = cosTheta * cosTheta;
   float SinTheta2 = 1.0 - CosTheta2;
   float Eta2 = eta * eta;
   float Etak2 = etak * etak;

   float t0 = Eta2 - Etak2 - SinTheta2;
   float a2plusb2 = sqrt(t0 * t0 + 4.0 * Eta2 * Etak2);
   float t1 = a2plusb2 + CosTheta2;
   float a = sqrt(0.5f * (a2plusb2 + t0));
   float t2 = 2 * a * cosTheta;
   float Rs = (t1 - t2) / (t1 + t2);

   float t3 = CosTheta2 * a2plusb2 + SinTheta2 * SinTheta2;
   float t4 = t2 * SinTheta2;   
   float Rp = Rs * (t3 - t4) / (t3 + t4);

   return 0.5 * (Rp + Rs);
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

/*
Return the specular component for Blinn-Phong shading with a fresnel coefficient with a complex eta (a metal).
*/
vec4 getSpecular_EtaComplex_Phong_Fresnel(vec4 C, float I, vec4 n, vec4 L, vec4 V, float eta_in, float eta_out, float etak_out, float exponent) {
     float eta = eta_out / eta_in;
     float etak = etak_out / eta_in;
     
     vec4 H = normalize(V + L);
     float F = fresnel_dieletric_conductor(eta, etak, dot(H, L));

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
     // vec4 specular = getSpecular_EtaReal_Phong_Fresnel(C, lightIntensity, n, L, V, 1.0, eta, shininess);
     vec4 specular = getSpecular_EtaComplex_Phong_Fresnel(C, lightIntensity, n, L, V, 1.0, eta, etak, shininess);
     
     vec4 finalColor = ambiant + diffuse + specular;
     fragColor = finalColor;
}
