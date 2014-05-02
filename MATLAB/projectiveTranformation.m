% #define WIDTH 640
% #define HEIGHT 480
% #define MAXDEPTH 10000
% #define HALFWIDTH WIDTH/2
% #define HALFHEIGHT HEIGHT/2
% #define HFOV 1.01447 //fXtoZ
% #define VFOV 0.789809 // fYtoZ
% #define COEFX WIDTH/HFOV
% #define COEFY HEIGHT/VFOV
% #define HFOVTANHALF tan(HFOV/2)
% #define VFOVTANHALF tan(VFOV/2)
% 
% //--------------------------------------------------------------
% static inline ofPoint g_worldToProjective(const ofPoint& p){
%     // this assumes a 640 x 480 resolution as per above defines
%     
%     // X_res / x_to_z * X_RW / z + X_res/2
%     // aProjective[i].X = (XnFloat)fCoeffX * aRealWorld[i].X / aRealWorld[i].Z + nHalfXres;
%     // aProjective[i].Y = nHalfYres - (XnFloat)fCoeffY * aRealWorld[i].Y / aRealWorld[i].Z;
%     // aProjective[i].Z = aRealWorld[i].Z;
%     
%     ofPoint projective;
% 	projective.x = COEFX * p.x / p.z + HALFWIDTH;
%     projective.y = HALFHEIGHT - COEFY * p.y / p.z;
%     projective.z = p.z;
% 	return projective;
% }
% 
% //--------------------------------------------------------------
% static inline ofPoint g_worldToProjective(const XnVector3D& p){
% 	return g_worldToProjective(toOf(p));
% }
% 
% //--------------------------------------------------------------
% static inline ofPoint g_projectiveToWorld(const ofPoint& p){
%     // this assumes a 640 x 480 resolution as per above defines
%     
%     // X_RW = (X_proj / X_res - 1/2) * Z * x_to_z
%     // XnDouble fNormalizedX = (aProjective[i].X / outputMode.nXRes - 0.5);
%     // aRealWorld[i].X = (XnFloat)(fNormalizedX * aProjective[i].Z * fXToZ);
%     // XnDouble fNormalizedY = (0.5 - aProjective[i].Y / outputMode.nYRes);
%     // aRealWorld[i].Y = (XnFloat)(fNormalizedY * aProjective[i].Z * fYToZ);
%     // aRealWorld[i].Z = aProjective[i].Z;
% 
%     ofPoint world;
%     world.x = (p.x / WIDTH - 0.5) * p.z * HFOV;
%     world.y = (0.5 - p.y / HEIGHT) * p.z * VFOV;
%     world.z = p.z;
% 	return world;
% }