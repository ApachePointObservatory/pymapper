/*=============================================================================
  Copyright (C) 2014 Allied Vision Technologies.  All Rights Reserved.

  Redistribution of this file, in original or modified form, without
  prior written consent of Allied Vision Technologies is prohibited.

-------------------------------------------------------------------------------

  File:        AsynchronousGrabWrite.c

  Description: The AsynchronousGrabWrite example will grab images asynchronously
               using VimbaC.

-------------------------------------------------------------------------------

  THIS SOFTWARE IS PROVIDED BY THE AUTHOR "AS IS" AND ANY EXPRESS OR IMPLIED
  WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF TITLE,
  NON-INFRINGEMENT, MERCHANTABILITY AND FITNESS FOR A PARTICULAR  PURPOSE ARE
  DISCLAIMED.  IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT,
  INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
  (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
  AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR
  TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
  OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

=============================================================================*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef WIN32
    #include <windows.h>
#else
    #include <unistd.h>
    #include <time.h>
    #include <pthread.h>
#endif

#include <VimbaC/Include/VimbaC.h>

#include "VmbTransform.h"

#include <AsynchronousGrabWrite.h>
#include <Bitmap.h>
#include "Common/PrintVimbaVersion.h"
#include "Common/DiscoverGigECameras.h"
#include "fitsio.h"


enum
{
    NUM_FRAMES  = 3
};

VmbBool_t       g_bVimbaStarted             = VmbBoolFalse;     // Remember if Vimba is started
VmbBool_t       g_bStreaming                = VmbBoolFalse;     // Remember if Vimba is streaming
VmbBool_t       g_bAcquiring                = VmbBoolFalse;     // Remember if Vimba is acquiring
VmbHandle_t     g_CameraHandle              = NULL;             // A handle to our camera
VmbFrame_t      g_Frames[NUM_FRAMES];                           // The frames we capture into
FrameInfos      g_eFrameInfos               = FrameInfos_Off;   // Remember if we should print out frame infos
VmbBool_t       g_bEnableColorProcessing    = VmbBoolFalse;     // enables color processing for frames
double          g_dFrameTime                = 0.0;              // Timestamp of last frame
VmbBool_t       g_bFrameTimeValid           = VmbBoolFalse;     // Remember if there was a last timestamp
VmbUint64_t     g_nFrameID                  = 0;                // ID of last frame
VmbBool_t       g_bFrameIDValid             = VmbBoolFalse;     // Remember if there was a last ID
#ifdef WIN32
double          g_dFrequency                = 0.0;              //Frequency of tick counter in Win32
#else
#endif

char fileNameBuf[30]; //will hold updated filenames (with numbers attached)
int frameNumber = 1; //iterator
double frameTimeBegin;
double lastFrameTime;

#ifdef WIN32
    HANDLE          g_hMutex = INVALID_HANDLE_VALUE;

    VmbBool_t CreateApiLock()
    {
        if( INVALID_HANDLE_VALUE != g_hMutex)
        {
            DestroyApiLock();
        }
        g_hMutex = CreateMutex( NULL, FALSE, NULL );
        return INVALID_HANDLE_VALUE != g_hMutex;
    }
    void DestroyApiLock()
    {
        if( INVALID_HANDLE_VALUE != g_hMutex)
        {
            CloseHandle( g_hMutex );
        }
    }
    VmbBool_t AquireApiLock()
    {
        if( WAIT_OBJECT_0 == WaitForSingleObject( g_hMutex, INFINITE ) )
        {
            return VmbBoolTrue;
        }
        return VmbBoolFalse;
    }
    void ReleaseApiLock()
    {
        ReleaseMutex( g_hMutex );
    }
#else
    pthread_mutex_t g_Mutex = PTHREAD_MUTEX_INITIALIZER;
    VmbBool_t CreateApiLock()
    {
        return VmbBoolTrue;
    }
    void DestroyApiLock()
    {
    }
    VmbBool_t AquireApiLock()
    {
        if(0 == pthread_mutex_lock( &g_Mutex ) )
        {
            return VmbBoolTrue;
        }
        return VmbBoolFalse;
    }
    void ReleaseApiLock()
    {
        pthread_mutex_unlock( &g_Mutex );
    }
#endif

// Method: ConfigureCamera
//
// Purpose: set parameters on the camera
VmbError_t  ConfigureCamera( VmbHandle_t pCameraHandle)
{
    VmbError_t          err                 = VmbErrorSuccess;      // The function result

    const char* PixelFormat = "Mono8";
    int Width = 960;
    int Height = 960;
    int OffsetX= 160; // offset columns such that the square image region is centered on the plate.
    int BinningHorizontal = 1;
    int BinningVertical = 1;
    int GVSPPacketSize = 500;
    double ExposureTimeAbs = 12000;
    const char* ExposureMode = "Timed";

    // Set ExposureTimeAbs
    err = VmbFeatureFloatSet( pCameraHandle, "ExposureTimeAbs", ExposureTimeAbs);
    if (VmbErrorSuccess == err)
    {
        printf("Set ExposureTimeAbs to %f\n", ExposureTimeAbs);
    }
    else
    {
        printf("Failed to set ExposureTimeAbs.  Err Code: %i\n", err);
    }

    // Set OffsetX
    err = VmbFeatureIntSet( pCameraHandle, "OffsetX", OffsetX);
    if (VmbErrorSuccess == err)
    {
        printf("Set OffsetX to %i\n", OffsetX);
    }
    else
    {
        printf("Failed to set OffsetX.  Err Code: %i\n", err);
    }

    // Set ExposureMode
    err = VmbFeatureEnumSet( pCameraHandle, "ExposureMode", ExposureMode);
    if (VmbErrorSuccess == err)
    {
        printf("Set ExposureMode to %s\n", ExposureMode);
    }
    else
    {
        printf("Failed to set ExposureTimeAbs.  Err Code: %i\n", err);
    }

    // Set GVSPPacketSize
    err = VmbFeatureIntSet( pCameraHandle, "GVSPPacketSize", GVSPPacketSize);
    if (VmbErrorSuccess == err)
    {
        printf("Set GVSPPacketSize to %i\n", GVSPPacketSize);
    }
    else
    {
        printf("Failed to set GVSPPacketSize.  Err Code: %i\n", err);
    }


    // Set pixel format
    err = VmbFeatureEnumSet( pCameraHandle, "PixelFormat", PixelFormat);
    if (VmbErrorSuccess == err)
    {
        printf("Set PixelFormat to %s\n", PixelFormat);
    }
    else
    {
        printf("Failed to set PixelFormat.  Err Code: %i\n", err);
    }

    // Set window width
    err = VmbFeatureIntSet( pCameraHandle, "Width", Width);
    if (VmbErrorSuccess == err)
    {
        printf("Set Width to %i\n", Width);
    }
    else
    {
        printf("Failed to set Width.  Err Code: %i\n", err);
    }

    // Set window width
    err = VmbFeatureIntSet( pCameraHandle, "Height", Height);
    if (VmbErrorSuccess == err)
    {
        printf("Set Height to %i\n", Height);
    }
    else
    {
        printf("Failed to set Height.  Err Code: %i\n", err);
    }

    // Set binning
    err = VmbFeatureIntSet( pCameraHandle, "BinningHorizontal", BinningHorizontal);
    if (VmbErrorSuccess == err)
    {
        printf("Set BinningHorizontal to %i\n", BinningHorizontal);
    }
    else
    {
        printf("Failed to set BinningHorizontal.  Err Code: %i\n", err);
    }
    err = VmbFeatureIntSet( pCameraHandle, "BinningVertical", BinningVertical);
    if (VmbErrorSuccess == err)
    {
        printf("Set BinningVertical to %i\n", BinningVertical);
    }
    else
    {
        printf("Failed to set BinningVertical.  Err Code: %i\n", err);
    }

    return err;

}

void printerror( int status)
{
    /*****************************************************/
    /* Print out cfitsio error messages and exit program */
    /*****************************************************/


    if (status)
    {
       fits_report_error(stderr, status); /* print error report */

       exit( status );    /* terminate the program, returning error status */
    }
    return;
}

    //
// Method: ProcessFrame
//
// Purpose: save image as bitmap
//
// Parameters:
// [in] pFrame frame to process data might be destroyed dependent on transform function used
//
// https://heasarc.gsfc.nasa.gov/docs/software/fitsio/quick/node7.html

//
// Method: GetTime
//
// Purpose: get time indicator
//
// Returns: time indicator in seconds for differential measurements
double GetTime()
{
#ifdef WIN32
    LARGE_INTEGER nCounter;
    QueryPerformanceCounter( &nCounter );
    return ( (double)nCounter.QuadPart ) / g_dFrequency;
#else
    struct timespec now;
    clock_gettime( CLOCK_REALTIME, &now );
    return ( (double)now.tv_sec ) + ( (double)now.tv_nsec ) / 1000000000.0;
#endif //WIN32
}


VmbError_t ProcessFrame( VmbFrame_t * pFrame )
{
    sprintf(fileNameBuf, "img%06d.fits", frameNumber);
    /******************************************************/
    /* Create a FITS primary array containing a 2-D image */
    /******************************************************/

    fitsfile *fptr;       /* pointer to the FITS file, defined in fitsio.h */
    int status, ii, jj;
    long  fpixel, nelements;
    double frameTime, fps;
    // unsigned short *array[960];

    /* initialize FITS image parameters */
    int bitpix   =  BYTE_IMG; //USHORT_IMG;
    long naxis    =   2;  /* 2-dimensional image                            */
    long naxes[2] = { 960, 960 };

    // remove(filename);               /* Delete old file if it already exists */

    status = 0;         /* initialize status before calling fitsio routines */

    // get timestamp stuff
    if (frameNumber == 1){
        frameTimeBegin = GetTime();
        lastFrameTime = 0;
        frameTime = 0;
        fps = 0;
    }
    else{

        frameTime = GetTime() - frameTimeBegin;
        fps =  1.0 / (frameTime - lastFrameTime);
    }
    if (fits_create_file(&fptr, fileNameBuf, &status)) /* create new FITS file */
         printerror( status );           /* call printerror if error occurs */

    /* write the required keywords for the primary array image.     */
    /* Since bitpix = USHORT_IMG, this will cause cfitsio to create */
    /* a FITS image with BITPIX = 16 (signed short integers) with   */
    /* BSCALE = 1.0 and BZERO = 32768.  This is the convention that */
    /* FITS uses to store unsigned integers.  Note that the BSCALE  */
    /* and BZERO keywords will be automatically written by cfitsio  */
    /* in this case.                                                */

    if ( fits_create_img(fptr,  bitpix, naxis, naxes, &status) )
         printerror( status );


    fpixel = 1;                               /* first pixel to write      */
    nelements = naxes[0] * naxes[1];          /* number of pixels to write */

    /* write the array of unsigned integers to the FITS file */
    if ( fits_write_img(fptr, TBYTE, fpixel, nelements, pFrame->buffer, &status) )
    // if ( fits_write_img(fptr, TUSHORT, fpixel, nelements, array[0], &status) )
        printerror( status );

    // free( array[0] );  /* free previously allocated memory */

    /* write another optional keyword to the header */
    /* Note that the ADDRESS of the value is passed in the routine */
    if ( fits_update_key(fptr, TDOUBLE, "TSTAMP", &frameTime,
         "timestamp of frame (seconds)", &status) )
         printerror( status );

    // printf("timestamp: %.4f   fps: %.4f", frameTime, fps);
    if ( fits_update_key(fptr, TDOUBLE, "FPS", &fps,
         "frames per second", &status) )
         printerror( status );

    if ( fits_update_key(fptr, TUINT, "FNum", &frameNumber,
         "frame number", &status) )
         printerror( status );

    if ( fits_close_file(fptr, &status) )                /* close the file */
         printerror( status );
    if (frameNumber % 100 == 0){
        printf( "fits successfully written to file \"%s\"\n", fileNameBuf );
    }
    frameNumber ++;
    lastFrameTime = frameTime;
    // return;
}



// VmbError_t ProcessFrame( VmbFrame_t * pFrame)
// {
//     AVTBitmap           bitmap;                                 // The bitmap we create


//     // check if we can get data
//     if( NULL == pFrame || NULL == pFrame->buffer )
//     {
//         printf("%s error invalid frame\n", __FUNCTION__);
//         return VmbErrorBadParameter;
//     }
//     bitmap.bufferSize = pFrame->imageSize;
//     bitmap.width = pFrame->width;
//     bitmap.height = pFrame->height;
//     bitmap.colorCode = ColorCodeMono8;
//     // Create the bitmap
//     if ( 0 == AVTCreateBitmap( &bitmap, pFrame->buffer ))
//     {
//         printf( "Could not create bitmap.\n" );
//     }
//     else
//     {
//         // Save the bitmap
//         // pFileName = "myFileName%i", ii;
//         sprintf(fileNameBuf, "img%i.bmp", frameNumber);
//         if ( 0 == AVTWriteBitmapToFile( &bitmap, fileNameBuf ))
//         {
//             printf( "Could not write bitmap to file.\n" );
//         }
//         else
//         {
//             printf( "Bitmap successfully written to file \"%s\"\n", fileNameBuf );
//             // Release the bitmap's buffer
//             frameNumber ++;
//             if ( 0 == AVTReleaseBitmap( &bitmap ))
//             {
//                 printf( "Could not release the bitmap.\n" );
//             }
//         }
//     }
// }


// //
// // Method: ProcessFrame
// //
// // Purpose: convert frames to RGB24 format and apply color processing if desired
// //
// // Parameters:
// // [in] pFrame frame to process data might be destroyed dependent on transform function used
// //
// VmbError_t ProcessFrame( VmbFrame_t * pFrame)
// {
//     VmbError_t          Result              = VmbErrorSuccess;  // result of function
//     VmbUint32_t         Width               = 0;                // will later hold the frame width
//     VmbUint32_t         Height              = 0;                // will later hold the frame height
//     VmbImage            SourceImage;                            // source image struct to pass to image transform
//     VmbImage            DestinationImage;                       // destination image struct to pass to image transform
//     VmbRGB8_t*          DestinationBuffer   = NULL;             // destination image buffer
//     VmbTransformInfo    TransformInfo;                          // if desired the transform information is constructed here
//     VmbUint32_t         TransformInfoCount  = 0;                // if color processing is desired this will be set
//     // check if we can get data
//     if( NULL == pFrame || NULL == pFrame->buffer )
//     {
//         printf("%s error invalid frame\n", __FUNCTION__);
//         return VmbErrorBadParameter;
//     }
//     // init local variables for frame width and height
//     Width   = pFrame->width;
//     Height  = pFrame->height;
//     if( g_bEnableColorProcessing == VmbBoolTrue )   // if color processing is desired set the transform matrix
//     {
//         const static VmbFloat_t matrix[9] = {   0.0, 0.0, 1.0,                  // matrix to swap red and blue component
//                                                 0.0, 1.0, 0.0,
//                                                 1.0, 0.0, 0.0 };
//         Result = VmbSetColorCorrectionMatrix3x3(matrix, &TransformInfo);        // initialize transform info
//         if( VmbErrorSuccess != Result)
//         {
//             printf("%s error could not set transform matrix; Error: %d\n", __FUNCTION__, Result);
//             return Result;
//         }
//         TransformInfoCount = 1;
//     }
//     // set the struct size to image
//     SourceImage.Size        = sizeof( SourceImage );
//     // set the image information from the frames pixel format and size
//     Result                  = VmbSetImageInfoFromPixelFormat( pFrame->pixelFormat, Width, Height, &SourceImage );
//     if( VmbErrorSuccess != Result)
//     {
//         printf( "%s error could not set source image info; Error: %d\n", __FUNCTION__, Result);
//         return Result;
//     }
//     // the frame buffer will be the images data buffer
//     SourceImage.Data = pFrame->buffer;
//     // set size for destination image
//     DestinationImage.Size   = sizeof( DestinationImage );
//     // set destination image info from frame size and string for RGB8 (rgb24)
//     Result                  = VmbSetImageInfoFromString( "RGB8", 4, Width, Height, &DestinationImage );
//     if( VmbErrorSuccess != Result)
//     {
//         printf("%s error could not set destination image info; Error: %d\n", __FUNCTION__, Result);
//         return Result;
//     }
//     // allocate buffer for destination image size is width * height * size of rgb pixel
//     DestinationBuffer       = (VmbRGB8_t*) malloc( Width*Height*sizeof( VmbRGB8_t) );
//     if( NULL == DestinationBuffer)
//     {
//         printf("%s error could not allocate rgb buffer for width: %d and height: %d\n", __FUNCTION__, Width, Height);
//         return VmbErrorResources;
//     }
//     // set the destination buffer to the data buffer of the image
//     DestinationImage.Data = DestinationBuffer;
//     // transform source to destination if color processing was enabled TransformInfoCount is 1 otherwise TransformInfo will be ignored
//     Result = VmbImageTransform( &SourceImage, &DestinationImage, &TransformInfo, TransformInfoCount );
//     // print first rgb pixel
//     printf("R: %d\tG: %d\tB: %d\n", DestinationBuffer->R, DestinationBuffer->G, DestinationBuffer->B);
//     // clean image buffer
//     free( DestinationBuffer );
//     return -Result;
// }


//
// Method: FrameCallback
//
// Purpose: called from Vimba if a frame is ready for user processing
//
// Parameters:
//
// [in] handle to camera that supplied the frame
// [in] pointer to frame structure that can hold valid data
//
void VMB_CALL FrameCallback( const VmbHandle_t cameraHandle, VmbFrame_t* pFrame )
{
    //
    // from here on the frame is under user control until returned to Vimba by re queuing it
    // if you want to have smooth streaming keep the time you hold the frame short
    //
    VmbBool_t       bShowFrameInfos     = VmbBoolFalse;         // showing frame infos
    double          dFPS                = 0.0;                  // frames per second calculated
    VmbBool_t       bFPSValid           = VmbBoolFalse;         // indicator if fps calculation was valid
    double          dFrameTime          = 0.0;                  // reference time for frames
    double          dTimeDiff           = 0.0;                  // time difference between frames
    VmbUint64_t     nFramesMissing      = 0;                    // number of missing frames

    // Ensure that a frame callback is not interrupted by a VmbFrameRevoke during shutdown
    AquireApiLock();

    if( FrameInfos_Off != g_eFrameInfos )
    {
        if( FrameInfos_Show == g_eFrameInfos )
        {
            bShowFrameInfos = VmbBoolTrue;
        }

        if( VmbFrameFlagsFrameID & pFrame->receiveFlags )
        {
            if( g_bFrameIDValid )
            {
                if( pFrame->frameID != ( g_nFrameID + 1 ) )
                {
                    // get difference between current frame and last received frame to calculate missing frames
                    nFramesMissing = pFrame->frameID - g_nFrameID - 1;
                    if( 1 == nFramesMissing )
                    {
                        printf("%s 1 missing frame detected\n", __FUNCTION__);
                    }
                    else
                    {
                        printf("%s error %llu missing frames detected\n",__FUNCTION__, nFramesMissing);
                    }
                }
            }
            g_nFrameID      = pFrame->frameID;          // store current frame id to calculate missing frames in the next calls
            g_bFrameIDValid = VmbBoolTrue;

            dFrameTime = GetTime();                     // get current time to calculate frames per second
            if(     ( g_bFrameTimeValid )               // only if the last time was valid
                &&  ( 0 == nFramesMissing ) )           // and the frame is not missing
            {
                dTimeDiff = dFrameTime - g_dFrameTime;  // build time difference with last frames time
                if( dTimeDiff > 0.0 )
                {
                    dFPS        = 1.0 / dTimeDiff;
                    bFPSValid   = VmbBoolTrue;
                }
                else
                {
                    bShowFrameInfos = VmbBoolTrue;
                }
            }
            // store time for fps calculation in the next call
            g_dFrameTime        = dFrameTime;
            g_bFrameTimeValid   = VmbBoolTrue;
        }
        else
        {
            bShowFrameInfos     = VmbBoolTrue;
            g_bFrameIDValid     = VmbBoolFalse;
            g_bFrameTimeValid   = VmbBoolFalse;
        }
        // test if the frame is complete
        if( VmbFrameStatusComplete != pFrame->receiveStatus )
        {
            bShowFrameInfos = VmbBoolTrue;
        }
    }

    if( bShowFrameInfos )
    {
        printf("Frame ID:");
        if( VmbFrameFlagsFrameID & pFrame->receiveFlags )
        {

            printf( "%llu", pFrame->frameID );
        }
        else
        {
            printf( "?" );
        }

        printf( " Status:" );
        switch( pFrame->receiveStatus )
        {
        case VmbFrameStatusComplete:
            printf( "Complete" );
            break;

        case VmbFrameStatusIncomplete:
            printf( "Incomplete" );
            break;

        case VmbFrameStatusTooSmall:
            printf( "Too small" );
            break;

        case VmbFrameStatusInvalid:
            printf( "Invalid" );
            break;

        default:
            printf( "?" );
            break;
        }

        printf( " Size:" );
        if( VmbFrameFlagsDimension & pFrame->receiveFlags )
        {
            printf( "%ux%u", pFrame->width, pFrame->height );
        }
        else
        {
            printf( "?x?" );
        }

        printf( " Format:0x%08X", pFrame->pixelFormat );

        printf( " FPS:" );
        if( bFPSValid )
        {
            printf( "%.2f", dFPS );
        }
        else
        {
            printf( "?" );
        }

        printf( "\n" );
    }
    // goto image processing
    ProcessFrame( pFrame );

    fflush( stdout );
    // requeue the frame so it can be filled again
    VmbCaptureFrameQueue( cameraHandle, pFrame, &FrameCallback );

    ReleaseApiLock();
}

//
// Method StartContinuousImageAcquisition
//
// Purpose: starts image acquisition on a given camera
//
// Parameters:
//
// [in]     pCameraId               zero terminated C string with the camera id for the camera to be used
// [in]     eFrameInfos             enumeration value for the frame infos to show for received frames
// [in]     bEnableColorProcessing  toggle for enabling image processing, in this case just swapping red with blue
//
// Note: Vimba has to be uninitialized and the camera has to allow access mode full
//
VmbError_t StartContinuousImageAcquisition( const char* pCameraID, FrameInfos eFrameInfos, VmbBool_t bEnableColorProcessing)
{
    VmbError_t          err                 = VmbErrorSuccess;      // The function result
    VmbCameraInfo_t     *pCameras           = NULL;                 // A list of camera details
    VmbUint32_t         nCount              = 0;                    // Number of found cameras
    VmbUint32_t         nFoundCount         = 0;                    // Change of found cameras
    VmbAccessMode_t     cameraAccessMode    = VmbAccessModeFull;    // We open the camera with full access
    VmbBool_t           bIsCommandDone      = VmbBoolFalse;         // Has a command finished execution
    VmbInt64_t          nPayloadSize        = 0;                    // The size of one frame
    int                 i                   = 0;                    // Counting variable
#ifdef WIN32
    LARGE_INTEGER       nFrequency;
#endif //WIN32

    if( !g_bVimbaStarted )
    {
        // initialize global state
        g_bStreaming                = VmbBoolFalse;
        g_bAcquiring                = VmbBoolFalse;
        g_CameraHandle              = NULL;
        memset( g_Frames, 0, sizeof( g_Frames ));
        g_dFrameTime                = 0.0;
        g_bFrameTimeValid           = VmbBoolFalse;
        g_nFrameID                  = 0;
        g_bFrameIDValid             = VmbBoolFalse;
        g_eFrameInfos               = eFrameInfos;
        g_bEnableColorProcessing    = bEnableColorProcessing;


#ifdef WIN32
        QueryPerformanceFrequency( &nFrequency );
        g_dFrequency        = (double)nFrequency.QuadPart;
#endif  //WIN32
        // Startup Vimba
        err = VmbStartup();
        // Print the version of Vimba
        PrintVimbaVersion();
        if ( VmbErrorSuccess == err )
        {
            g_bVimbaStarted = VmbBoolTrue;

            // Is Vimba connected to a GigE transport layer?
            DiscoverGigECameras();

            // If no camera ID was provided use the first camera found
            if ( NULL == pCameraID )
            {
                // Get the amount of known cameras
                err = VmbCamerasList( NULL, 0, &nCount, sizeof *pCameras );
                if (    VmbErrorSuccess == err
                     && 0 != nCount )
                {
                    pCameras = (VmbCameraInfo_t*)malloc( nCount * sizeof( *pCameras ));
                    if ( NULL != pCameras )
                    {
                        // Actually query all static details of all known cameras without having to open the cameras
                        // If a new camera was connected since we queried the amount of cameras (nFoundCount > nCount) we can ignore that one
                        err = VmbCamerasList( pCameras, nCount, &nFoundCount, sizeof *pCameras );
                        if (    VmbErrorSuccess != err
                             && VmbErrorMoreData != err )
                        {
                            printf( "%s Could not list cameras. Error code: %d\n", __FUNCTION__, err );
                        }
                        else
                        {
                            // Use the first camera
                            if( nFoundCount != 0)
                            {
                                pCameraID = pCameras[0].cameraIdString;
                            }
                            else
                            {
                                err = VmbErrorNotFound;
                                printf( "%s camera lost. Error code: %d\n", __FUNCTION__, err );
                                pCameraID = NULL;
                            }
                        }

                        free( pCameras );
                        pCameras = NULL;
                    }
                    else
                    {
                        printf( "%s Could not allocate camera list.\n", __FUNCTION__ );
                    }
                }
                else
                {
                    printf( "%s Could not list cameras or no cameras present. Error code: %d\n", __FUNCTION__, err );
                }
            }

            if ( NULL != pCameraID )
            {
                // Open camera
                err = VmbCameraOpen( pCameraID, cameraAccessMode, &g_CameraHandle );
                if ( VmbErrorSuccess == err )
                {
                    printf("Opening camera with ID: %s\n", pCameraID);
                    // (In this example we do not test whether this cam actually is a GigE cam)
                    if ( VmbErrorSuccess == VmbFeatureCommandRun( g_CameraHandle, "GVSPAdjustPacketSize" ))
                    {
                        do
                        {
                            if ( VmbErrorSuccess != VmbFeatureCommandIsDone(    g_CameraHandle,
                                                                                "GVSPAdjustPacketSize",
                                                                                &bIsCommandDone ))
                            {
                                break;
                            }
                        } while ( VmbBoolFalse == bIsCommandDone );
                    }

                    if ( VmbErrorSuccess == err )
                    {
                        err = ConfigureCamera( g_CameraHandle );
                        printf("ConfiguredCamera err = %i\n", err);
                        // Evaluate frame size
                        err = VmbFeatureIntGet( g_CameraHandle, "PayloadSize", &nPayloadSize );
                        if ( VmbErrorSuccess == err )
                        {
                            for(i = 0; i < NUM_FRAMES; i++)
                            {
                                g_Frames[i].buffer        = (unsigned char*)malloc( (VmbUint32_t)nPayloadSize );
                                if( NULL == g_Frames[i].buffer )
                                {
                                    err = VmbErrorResources;
                                    break;
                                }
                                g_Frames[i].bufferSize = (VmbUint32_t)nPayloadSize;

                                // Announce Frame
                                err = VmbFrameAnnounce( g_CameraHandle, &g_Frames[i], (VmbUint32_t)sizeof( VmbFrame_t ));
                                if ( VmbErrorSuccess != err )
                                {
                                    free( g_Frames[i].buffer );
                                    memset( &g_Frames[i], 0, sizeof( VmbFrame_t ));
                                    break;
                                }
                            }

                            if ( VmbErrorSuccess == err )
                            {
                                // Start Capture Engine
                                err = VmbCaptureStart( g_CameraHandle );
                                if ( VmbErrorSuccess == err )
                                {
                                    g_bStreaming = VmbBoolTrue;
                                    for( i = 0; i < NUM_FRAMES; i++ )
                                    {
                                        // Queue Frame
                                        err = VmbCaptureFrameQueue( g_CameraHandle, &g_Frames[i], &FrameCallback );
                                        if ( VmbErrorSuccess != err )
                                        {
                                            break;
                                        }
                                    }

                                    if ( VmbErrorSuccess == err )
                                    {
                                        // Start Acquisition
                                        err = VmbFeatureCommandRun( g_CameraHandle,"AcquisitionStart" );
                                        if ( VmbErrorSuccess == err )
                                        {
                                            g_bAcquiring = VmbBoolTrue;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }

            if( VmbErrorSuccess != err )
            {
                StopContinuousImageAcquisition();
            }
        }
    }
    else
    {
        err = VmbErrorOther;
    }

    return err;
}

//
// Method: StopContinuousImageAcquisition
//
// Purpose: stops image acquisition that was started with StartContinuousImageAcquisition
//
void StopContinuousImageAcquisition()
{
    int i = 0;

    if( g_bVimbaStarted )
    {
        if( NULL != g_CameraHandle )
        {
            if( g_bAcquiring )
            {
                // Stop Acquisition
                VmbFeatureCommandRun( g_CameraHandle, "AcquisitionStop" );
                g_bAcquiring = VmbBoolFalse;
            }

            if( g_bStreaming )
            {
                // Stop Capture Engine
                VmbCaptureEnd( g_CameraHandle );
                g_bStreaming = VmbBoolFalse;
            }

            // Flush the capture queue
            VmbCaptureQueueFlush( g_CameraHandle );

            // Ensure that revoking is not interrupted by a dangling frame callback
            AquireApiLock();
            for( i = 0; i < NUM_FRAMES; i++ )
            {
                if( NULL != g_Frames[i].buffer )
                {
                    VmbFrameRevoke( g_CameraHandle, &g_Frames[i] );
                    free( g_Frames[i].buffer );
                    memset( &g_Frames[i], 0, sizeof( VmbFrame_t ));
                }
            }
            ReleaseApiLock();
            // Close camera
            VmbCameraClose ( g_CameraHandle );
            g_CameraHandle = NULL;
        }
        VmbShutdown();
        g_bVimbaStarted = VmbBoolFalse;
#ifdef WIN32
        CloseHandle( g_hMutex );
#else
        pthread_mutex_destroy(&g_Mutex);
#endif
    }
}

