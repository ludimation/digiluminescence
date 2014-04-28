#include "cinder/app/AppBasic.h"
#include "cinder/gl/gl.h"
#include "cinder/gl/Texture.h"

#include "CinderFreenect.h"

using namespace ci;
using namespace ci::app;
using namespace std;

class digiluminescenceApp : public AppBasic {
  public:
	void prepareSettings( Settings* settings );
	void setup();
	void update();
	void draw();
	
	KinectRef		mKinect;
	gl::Texture		mColorTexture, mDepthTexture;	
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

CINDER_APP_BASIC( digiluminescenceApp, RendererGl )
