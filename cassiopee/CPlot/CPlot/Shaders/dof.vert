#version 130
// Depth of field shader
void main()
{
   gl_TexCoord[0] = gl_MultiTexCoord0;
   gl_Position = ftransform();
}
