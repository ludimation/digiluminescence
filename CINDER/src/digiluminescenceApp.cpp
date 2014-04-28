// From KinectBasic.cpp
#include "cinder/app/AppBasic.h"
#include "cinder/gl/gl.h"
#include "cinder/gl/Texture.h"

#include "CinderFreenect.h"

using namespace ci;
using namespace ci::app;
using namespace std;

// From SkeletonApp.cpp
#include "cinder/Camera.h"
#include "cinder/ImageIo.h"
#include "cinder/params/Params.h"
#include "cinder/Utilities.h"

class digiluminescenceApp : public AppBasic {
public:
	void prepareSettings( Settings* settings );
	void setup();
	void update();
	void draw();
    
	gl::Texture		mColorTexture, mDepthTexture;
    
    // From SkeletonApp.cpp
    void	keyDown( ci::app::KeyEvent event );
	void	shutdown();

private:
    
    // From SkeletonApp.cpp
	// Kinect
	KinectRef		mKinect;
	uint32_t							mCallbackId;
    // need to try and find freenect/open NI equivalents to following lines if possible
    // -- Looks negative -- http://openkinect.org/wiki/FAQ#Does_libfreenect_have_any_skeleton_tracking_feature.3F
    //	std::vector<Skeleton>	mSkeletons;
    //	void								onSkeletonData( std::vector<Skeleton> skeletons,
    //                                                       const DeviceOptions &deviceOptions );
    
	// Camera
	ci::CameraPersp						mCamera;
    
	// Save screenshot
	void								screenShot();

    // TODO: find an addon that does skeleton trancking stuff (might need to use Open NI, but not sure that is being supported by Cinder at this point)
    
};

void digiluminescenceApp::prepareSettings( Settings* settings )
{
	settings->setWindowSize( 1280, 480 );
}

void digiluminescenceApp::setup()
{
	mKinect = Kinect::create();
}

void digiluminescenceApp::update()
{	
	if( mKinect->checkNewDepthFrame() )
		mDepthTexture = mKinect->getDepthImage();
	
	if( mKinect->checkNewVideoFrame() )
		mColorTexture = mKinect->getVideoImage();
}

void digiluminescenceApp::draw()
{
	gl::clear(); 
	gl::setMatricesWindow( getWindowWidth(), getWindowHeight() );
	if( mDepthTexture )
		gl::draw( mDepthTexture );
	if( mColorTexture )
		gl::draw( mColorTexture, Vec2i( 640, 0 ) );
}

// Handles key press
void digiluminescenceApp::keyDown( KeyEvent event )
{
    
	// Quit, toggle fullscreen
	switch ( event.getCode() ) {
        case KeyEvent::KEY_q:
            quit();
            break;
        case KeyEvent::KEY_f:
            setFullScreen( !isFullScreen() );
            break;
        case KeyEvent::KEY_SPACE:
           // screenShot();
            break;
	}
    
}

CINDER_APP_BASIC( digiluminescenceApp, RendererGl )
