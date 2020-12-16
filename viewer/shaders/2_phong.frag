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

vec4 getAmbiant(float k, vec4 C, float I) {
     return k * C * I;
}

vec4 getDiffuse(float k, vec4 C, float I, vec4 n, vec4 L) {
     return k * C * max(dot(n, L), 0.0) * I;
}

void main( void )
{
     vec4 C = vertColor;

     vec4 n = normalize(vertNormal);
     vec4 L = normalize(lightVector);

     // Ambiant 
     float kA = 0.2;
     vec4 ambiant = getAmbiant(kA, C, lightIntensity);
     
     // Diffuse
     float kD = 0.8;
     vec4 diffuse = getDiffuse(kD, C, lightIntensity, n, L);
     
     // Specular
     vec4 specular = vec4(0);
     
     vec4 finalColor = ambiant + diffuse + specular;
     fragColor = finalColor;
}
