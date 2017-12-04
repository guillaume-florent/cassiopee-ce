// Depth of field shader
uniform sampler2D FrameBuffer;
uniform sampler2D depthMap;

uniform int windowWidth, windowHeight;
uniform float focalDepth; // position de la focale
uniform float radius; // taille du rayon de blur

vec2 poisson0, poisson1, poisson2, poisson3, poisson4;
vec2 poisson5, poisson6, poisson7;
vec2 maxCoC = vec2(radius, 2.*radius);

vec2 pixelSizeHigh;

vec4 dof(vec2 coords)
{
     // pas de blur au point focal
     // puis lineaire avant et apres sur dist en depth
     float dist = 0.02/(2.*radius);
     vec4 finalColor;
     float discRadius, discRadiusLow, centerDepth;

     centerDepth = texture2D(depthMap, coords).r;
     centerDepth = (centerDepth - focalDepth)/dist;
     centerDepth = clamp(centerDepth,-1.,1.);
     centerDepth = centerDepth*0.5 + 0.5;

     discRadius = abs(centerDepth*maxCoC.y - maxCoC.x);
     discRadiusLow = discRadius * 0.5;

     finalColor = vec4(0.0,0.0,0.0,0.0);
     float finalBlur = 0.;
     float blur, tapBlur;
    
     // -- 
     vec2 coordHigh = coords + (pixelSizeHigh*poisson0*discRadius);
     vec4 tapHigh = texture2D(FrameBuffer, coordHigh);
     float blurHigh = texture2D(depthMap, coordHigh).r;
     blurHigh = (blurHigh - focalDepth)/dist;
     blurHigh = clamp(blurHigh,-1.,1.);
     blurHigh = blurHigh*0.5 + 0.5;
     tapBlur = abs(blurHigh*2.0 - 1.0);
     blur = ((blurHigh >= centerDepth) ? 1.0 : tapBlur);
     finalColor += tapHigh * blur;
     finalBlur += blur;
      
     // -- 
     coordHigh = coords + (pixelSizeHigh*poisson1*discRadius);
     tapHigh = texture2D(FrameBuffer, coordHigh);
     blurHigh = texture2D(depthMap, coordHigh).r;
     blurHigh = (blurHigh - focalDepth)/dist;
     blurHigh = clamp(blurHigh,-1.,1.);
     blurHigh = blurHigh*0.5 + 0.5;
     tapBlur = abs(blurHigh*2.0 - 1.0);
     blur = ((blurHigh >= centerDepth) ? 1.0 : tapBlur);
     finalColor += tapHigh * blur;
     finalBlur += blur;

     // -- 
     coordHigh = coords + (pixelSizeHigh*poisson2*discRadius);
     tapHigh = texture2D(FrameBuffer, coordHigh);
     blurHigh = texture2D(depthMap, coordHigh).r;
     blurHigh = (blurHigh - focalDepth)/dist;
     blurHigh = clamp(blurHigh,-1.,1.);
     blurHigh = blurHigh*0.5 + 0.5;
     tapBlur = abs(blurHigh*2.0 - 1.0);
     blur = ((blurHigh >= centerDepth) ? 1.0 : tapBlur);
     finalColor += tapHigh * blur;
     finalBlur += blur;

     // -- 
     coordHigh = coords + (pixelSizeHigh*poisson3*discRadius);
     tapHigh = texture2D(FrameBuffer, coordHigh);
     blurHigh = texture2D(depthMap, coordHigh).r;
     blurHigh = (blurHigh - focalDepth)/dist;
     blurHigh = clamp(blurHigh,-1.,1.);
     blurHigh = blurHigh*0.5 + 0.5;
     tapBlur = abs(blurHigh*2.0 - 1.0);
     blur = ((blurHigh >= centerDepth) ? 1.0 : tapBlur);
     finalColor += tapHigh * blur;
     finalBlur += blur;

     // -- 
     coordHigh = coords + (pixelSizeHigh*poisson4*discRadius);
     tapHigh = texture2D(FrameBuffer, coordHigh);
     blurHigh = texture2D(depthMap, coordHigh).r;
     blurHigh = (blurHigh - focalDepth)/dist;
     blurHigh = clamp(blurHigh,-1.,1.);
     blurHigh = blurHigh*0.5 + 0.5;
     tapBlur = abs(blurHigh*2.0 - 1.0);
     blur = ((blurHigh >= centerDepth) ? 1.0 : tapBlur);
     finalColor += tapHigh * blur;
     finalBlur += blur;

     // -- 
     coordHigh = coords + (pixelSizeHigh*poisson5*discRadius);
     tapHigh = texture2D(FrameBuffer, coordHigh);
     blurHigh = texture2D(depthMap, coordHigh).r;
     blurHigh = (blurHigh - focalDepth)/dist;
     blurHigh = clamp(blurHigh,-1.,1.);
     blurHigh = blurHigh*0.5 + 0.5;
     tapBlur = abs(blurHigh*2.0 - 1.0);
     blur = ((blurHigh >= centerDepth) ? 1.0 : tapBlur);
     finalColor += tapHigh * blur;
     finalBlur += blur;

     // -- 
     coordHigh = coords + (pixelSizeHigh*poisson6*discRadius);
     tapHigh = texture2D(FrameBuffer, coordHigh);
     blurHigh = texture2D(depthMap, coordHigh).r;
     blurHigh = (blurHigh - focalDepth)/dist;
     blurHigh = clamp(blurHigh,-1.,1.);
     blurHigh = blurHigh*0.5 + 0.5;
     tapBlur = abs(blurHigh*2.0 - 1.0);
     blur = ((blurHigh >= centerDepth) ? 1.0 : tapBlur);
     finalColor += tapHigh * blur;
     finalBlur += blur;

     // -- 
     coordHigh = coords + (pixelSizeHigh*poisson7*discRadius);
     tapHigh = texture2D(FrameBuffer, coordHigh);
     blurHigh = texture2D(depthMap, coordHigh).r;
     blurHigh = (blurHigh - focalDepth)/dist;
     blurHigh = clamp(blurHigh,-1.,1.);
     blurHigh = blurHigh*0.5 + 0.5;
     tapBlur = abs(blurHigh*2.0 - 1.0);
     blur = ((blurHigh >= centerDepth) ? 1.0 : tapBlur);
     finalColor += tapHigh * blur;
     finalBlur += blur;

     return finalColor/finalBlur;
}

void main(){

	pixelSizeHigh[0] = 1.0/float(windowWidth);
	pixelSizeHigh[1] = 1.0/float(windowHeight);

	// poisson-distributed positions
        poisson0 = vec2( 0.0, 0.0);
        poisson1 = vec2( 0.527837,-0.08586);
        poisson2 = vec2(-0.040088, 0.536087);
        poisson3 = vec2(-0.670445,-0.179949);
        poisson4 = vec2(-0.419418,-0.616039);
        poisson5 = vec2( 0.440453,-0.639399);
        poisson6 = vec2(-0.757088, 0.349334);
        poisson7 = vec2( 0.574619, 0.685879);

	gl_FragColor = dof(gl_TexCoord[0].xy);
}
