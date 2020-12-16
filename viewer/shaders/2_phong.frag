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

void main( void )
{
     vec4 C = vertColor;

     // Ambiant 
     float kA = 1;
     vec4 ambiant = getAmbiant(kA, C, lightIntensity);
     
     // Diffuse
     vec4 diffuse = vec4(0);
     
     // Specular
     vec4 specular = vec4(0);
     
     vec4 finalColor = ambiant + diffuse + specular;
     fragColor = finalColor;
}
