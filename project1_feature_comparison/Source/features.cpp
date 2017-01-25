#include <assert.h>
#include <math.h>
#include <FL/Fl.H>
#include <FL/Fl_Image.H>
#include "features.h"
#include "ImageLib/FileIO.h"

#define PI 3.14159265358979323846
using namespace std;

// Compute features of an image.
bool computeFeatures(CFloatImage &image, FeatureSet &features, int featureType) {
    // TODO: Instead of calling dummyComputeFeatures, write your own
    // feature computation routines and call them here.
    switch (featureType) {
        case 1:
            dummyComputeFeatures(image, features);
            break;
        case 2:
            ComputeHarrisFeatures(image, features);
            ComputeMOPSDescriptors(image, features);
            //ComputeWinDescriptors(image, features);
            //ComputeHoGDescriptors(image, features);
            break;
        default:
            return false;
    }
    
    // This is just to make sure the IDs are assigned in order, because
    // the ID gets used to index into the feature array.
    for (unsigned int i=0; i<features.size(); i++) {
        features[i].id = i+1;
    }
    
    return true;
}

// Perform a query on the database.  This simply runs matchFeatures on
// each image in the database, and returns the feature set of the best
// matching image.
bool performQuery(const FeatureSet &f, const ImageDatabase &db, int &bestIndex, vector<FeatureMatch> &bestMatches, double &bestScore, int matchType) {
    // Here's a nice low number.
    bestScore = -1e100;
    
    vector<FeatureMatch> tempMatches;
    double tempScore;
    
    for (unsigned int i=0; i<db.size(); i++) {
        if (!matchFeatures(f, db[i].features, tempMatches, tempScore, matchType)) {
            return false;
        }
        
        if (tempScore > bestScore) {
            bestIndex = i;
            bestScore = tempScore;
            bestMatches = tempMatches;
        }
    }
    
    return true;
}

// Match one feature set with another.
bool matchFeatures(const FeatureSet &f1, const FeatureSet &f2, vector<FeatureMatch> &matches, double &totalScore, int matchType) {
    // TODO: We have given you the ssd matching function, you must write your own
    // feature matching function for the ratio test.
    
    printf("\nMatching features.......\n");
    
    switch (matchType) {
        case 1:
            ssdMatchFeatures(f1, f2, matches, totalScore);
            return true;
        case 2:
            ratioMatchFeatures(f1, f2, matches, totalScore);
            return true;
        default:
            return false;
    }
}

// Evaluate a match using a ground truth homography.  This computes the
// average SSD distance between the matched feature points and
// the actual transformed positions.
double evaluateMatch(const FeatureSet &f1, const FeatureSet &f2, const vector<FeatureMatch> &matches, double h[9]) {
    double d = 0;
    int n = 0;
    
    double xNew;
    double yNew;
    
    unsigned int num_matches = matches.size();
    for (unsigned int i=0; i<num_matches; i++) {
        int id1 = matches[i].id1;
        int id2 = matches[i].id2;
        applyHomography(f1[id1-1].x, f1[id1-1].y, xNew, yNew, h);
        d += sqrt(pow(xNew-f2[id2-1].x,2)+pow(yNew-f2[id2-1].y,2));
        n++;
    }
    
    return d / n;
}

void addRocData(const FeatureSet &f1, const FeatureSet &f2, const vector<FeatureMatch> &matches, double h[9],vector<bool> &isMatch,double threshold,double &maxD) {
    double d = 0;
    
    double xNew;
    double yNew;
    
    unsigned int num_matches = matches.size();
    for (unsigned int i=0; i<num_matches; i++) {
        int id1 = matches[i].id1;
        int id2 = matches[i].id2;
        applyHomography(f1[id1-1].x, f1[id1-1].y, xNew, yNew, h);
        
        // Ignore unmatched points.  There might be a better way to
        // handle this.
        d = sqrt(pow(xNew-f2[id2-1].x,2)+pow(yNew-f2[id2-1].y,2));
        if (d<=threshold)
        {
            isMatch.push_back(1);
        }
        else
        {
            isMatch.push_back(0);
        }
        
        if (matches[i].score>maxD)
            maxD=matches[i].score;
    }
}

vector<ROCPoint> computeRocCurve(vector<FeatureMatch> &matches,vector<bool> &isMatch,vector<double> &thresholds)
{
    vector<ROCPoint> dataPoints;
    
    for (int i=0; i < (int)thresholds.size();i++)
    {
        //printf("Checking threshold: %lf.\r\n",thresholds[i]);
        int tp=0;
        int actualCorrect=0;
        int fp=0;
        int actualError=0;
        int total=0;
        
        int num_matches = (int) matches.size();
        for (int j=0;j < num_matches;j++)
        {
            if (isMatch[j])
            {
                actualCorrect++;
                if (matches[j].score<thresholds[i])
                {
                    tp++;
                }
            }
            else
            {
                actualError++;
                if (matches[j].score<thresholds[i])
                {
                    fp++;
                }
            }
            
            total++;
        }
        
        ROCPoint newPoint;
        //printf("newPoints: %lf,%lf",newPoint.trueRate,newPoint.falseRate);
        newPoint.trueRate=(double(tp)/actualCorrect);
        newPoint.falseRate=(double(fp)/actualError);
        //printf("newPoints: %lf,%lf",newPoint.trueRate,newPoint.falseRate);
        
        dataPoints.push_back(newPoint);
    }
    
    return dataPoints;
}



// Compute silly example features.  This doesn't do anything
// meaningful.
void dummyComputeFeatures(CFloatImage &image, FeatureSet &features) {
    CShape sh = image.Shape();
    Feature f;
    cout << "Computing dummy..." << endl;
    
    for (int y=0; y<sh.height; y++) {
        for (int x=0; x<sh.width; x++) {
            double r = image.Pixel(x,y,0);
            double g = image.Pixel(x,y,1);
            double b = image.Pixel(x,y,2);
            
            if ((int)(255*(r+g+b)+0.5) % 100  == 1) {
                // If the pixel satisfies this meaningless criterion,
                // make it a feature.
                
                f.type = 1;
                f.id += 1;
                f.x = x;
                f.y = y;
                
                f.data.resize(1);
                f.data[0] = r + g + b;
                
                features.push_back(f);
            }
        }
    }
}

void ComputeHarrisFeatures(CFloatImage &image, FeatureSet &features)
{
    
    cout << "Computing Harris:" << endl;
    //Create grayscale image used for Harris detection
    CFloatImage grayImage=ConvertToGray(image);
    
    CFloatImage Ix(grayImage.Shape().width,grayImage.Shape().height,1);
    CFloatImage Iy(grayImage.Shape().width,grayImage.Shape().height,1);
    
    CFloatImage angleImage(grayImage.Shape().width,grayImage.Shape().height,1);
    CFloatImage gradientImage(grayImage.Shape().width,grayImage.Shape().height,1);
    //Create image to store Harris values
    CFloatImage harrisImage(grayImage.Shape().width,grayImage.Shape().height,1);
    
    //Create image to store local maximum harris values as 1, other pixels 0
    CByteImage harrisMaxImage(grayImage.Shape().width,grayImage.Shape().height,1);
    
    //const int& h = grayImage.Shape().height;
    //const int& w = grayImage.Shape().width;
    
    
    //compute Harris values puts harris values at each pixel position in harrisImage.
    //You'll need to implement this function.
    computeHarrisValues(grayImage, harrisImage,Ix,Iy);
    
    for(int y=0;y<harrisMaxImage.Shape().height;y++){
        for(int x=0;x<harrisMaxImage.Shape().width;x++){
            angleImage.Pixel(x, y, 0)= atan2(Iy.Pixel(x, y, 0),Ix.Pixel(x, y, 0));
            gradientImage.Pixel(x, y, 0)=sqrt( Ix.Pixel(x, y, 0)*Ix.Pixel(x, y, 0)+Iy.Pixel(x, y, 0)*Iy.Pixel(x, y, 0) );
        }
    }
    
    // Threshold the harris image and compute local maxima.  You'll need to implement this function.
    int r_value=1;
    while(computeLocalMaxima(harrisImage,harrisMaxImage,r_value)>2000){
        r_value++;
    }
    // Prints out the harris image for debugging purposes
        CByteImage tmp(harrisImage.Shape());
        convertToByteImage(harrisImage, tmp);
        WriteFile(tmp, "/Users/heyu/Documents/cplusplus/PackageProject_1/harris2.tga");
        //WriteFile(harrisMaxImage, "harrisMax.tga");
    
    // TO DO--------------------------------------------------------------------
    //Loop through feature points in harrisMaxImage and create feature descriptor
    //for each point above a threshold
    
    
    
    
    for (int y=8;y<harrisMaxImage.Shape().height-8;y++) {
        for (int x=8;x<harrisMaxImage.Shape().width-8;x++) {
            
            // Skip over non-maxima
            if (harrisMaxImage.Pixel(x, y, 0) == 0)continue;//use local maximum only
            // Fill in feature with descriptor data here.
            Feature f;
            f.type = 2;
            f.id += 1;
            f.x = x;
            f.y = y;
            f.angleRadians=angleImage.Pixel(x, y, 0);
           
            features.push_back(f);
        }
    }
    //cout<<sizeof(features)<<endl;
}

void ComputeWinDescriptors(CFloatImage &image, FeatureSet &features){
    cout << "Computing Window..." << std::endl;
    // This image represents the window around the feature you need to compute to store as the feature descriptor
    const int windowSize = 5;
    CFloatImage destImage(windowSize, windowSize, 1);
    
    // Gray image
    CFloatImage grayImage = ConvertToGray(image);
    
    // Blur Image
    CFloatImage blurImage(grayImage.Shape().width, grayImage.Shape().height, 1);
    
    // Blur 5x5 Gaussian kernel
    int blurWindowSize = 5;
    CFloatImage blurKernel(blurWindowSize, blurWindowSize, 1);
    for(int i = 0; i < blurWindowSize; i++){
        for(int j = 0; j < blurWindowSize; j++){
            blurKernel.Pixel(j,i,0) = gaussian5x5[blurWindowSize * i + j];
        }
    }
    blurKernel.origin[0] = blurWindowSize/2; // x offset to the center of the kernel
    blurKernel.origin[1] = blurWindowSize/2; // y offset to the center of the kernel
    Convolve(grayImage, blurImage, blurKernel);
    
    for (vector<Feature>::iterator i = features.begin(); i != features.end(); i++) {
        Feature &f = *i;
        f.data.resize(windowSize * windowSize);
        for (int row = -windowSize/2; row <= +windowSize/2; row++ ){
            for (int col = -windowSize/2; col <= +windowSize/2; col++){
                f.data[(row+windowSize/2) * windowSize + col+windowSize/2] = grayImage.Pixel(f.x+row, f.y+col, 0);
            }
        }
    }
    
    
}


void ComputeHoGDescriptors(CFloatImage &image, FeatureSet &features)
{
    
    cout << "Computing HoG..." << std::endl;
    // This image represents the window around the feature you need to compute to store as the feature descriptor
    const int windowSize = 5;
    CFloatImage destImage(windowSize, windowSize, 1);
    
    // Gray image
    CFloatImage grayImage = ConvertToGray(image);
    CFloatImage Ix(grayImage.Shape().width,grayImage.Shape().height,1);
    CFloatImage Iy(grayImage.Shape().width,grayImage.Shape().height,1);
    CFloatImage gradientImage(grayImage.Shape().width,grayImage.Shape().height,1);
    CFloatImage angleImage(grayImage.Shape().width,grayImage.Shape().height,1);
    
    
    float Bx[3][3]={{-1,0,1},{-2,0,2},{-1,0,1}};
    float By[3][3]={{-1,-2,-1},{0,0,0},{1,2,1}};
    
    for (int x = 0; x < grayImage.Shape().width; x++) {
        for (int y = 0; y < grayImage.Shape().height; y++) {
            if(x>1&&y>1&&x+1<grayImage.Shape().width-1&&y+1<grayImage.Shape().height-1){
                for (int kx=0;kx<3;kx++){
                    for (int ky=0;ky<3;ky++){
                        Ix.Pixel(x,y,0)=Ix.Pixel(x,y,0)+grayImage.Pixel(x +kx-1,y +ky-1,0)*Bx[3-1-kx][3-1-ky];
                        Iy.Pixel(x,y,0)=Iy.Pixel(x,y,0)+grayImage.Pixel(x+kx-1,y+ky-1,0)*By[3-1-kx][3-1-ky];
                    }
                }
            }
        }
    }
    
    // Blur Image
    CFloatImage blurImage(image.Shape().width, image.Shape().height, 1);
    
    // Blur 5x5 Gaussian kernel
    int blurWindowSize = 5;
    CFloatImage blurKernel(blurWindowSize, blurWindowSize, 1);
    for(int i = 0; i < blurWindowSize; i++){
        for(int j = 0; j < blurWindowSize; j++){
            blurKernel.Pixel(j,i,0) = gaussian5x5[blurWindowSize * i + j];
        }
    }
    blurKernel.origin[0] = blurWindowSize/2; // x offset to the center of the kernel
    blurKernel.origin[1] = blurWindowSize/2; // y offset to the center of the kernel
    Convolve(grayImage, blurImage, blurKernel);
    
    for(int y=0;y<blurImage.Shape().height;y++){
        for(int x=0;x<blurImage.Shape().width;x++){
            angleImage.Pixel(x, y, 0)= atan2(Iy.Pixel(x, y, 0),Ix.Pixel(x, y, 0));
            gradientImage.Pixel(x, y, 0)=sqrt( Ix.Pixel(x, y, 0)*Ix.Pixel(x, y, 0)+Iy.Pixel(x, y, 0)*Iy.Pixel(x, y, 0) );
        }
    }
    
    for (vector<Feature>::iterator i = features.begin(); i != features.end(); i++) {
        Feature &f = *i;
        f.data.resize(windowSize * windowSize);
        
        
        vector<float> value(25,0);
        vector<float> direction(32,0);
        int windowsize=8;
        int startx[4]={f.x-windowsize/2,f.x,f.x-windowsize/2,f.x};
        int starty[4]={f.y-windowsize/2,f.y-windowsize/2,f.y,f.y};
        float total_len=0.0;
        for(int cell=0;cell<4;cell++){
            
            for(int i=startx[cell];i<startx[cell]+windowsize/2;i++){
                for(int j=starty[cell];j<starty[cell]+windowsize/2;j++){
                    
                    angleImage.Pixel(f.x+i,f.y+j,0)=angleImage.Pixel(f.x+i,f.y+j,0)-angleImage.Pixel(f.x,f.y,0);
                    if (angleImage.Pixel(f.x+i,f.y+j,0)>=-PI&&angleImage.Pixel(f.x+i,f.y+j,0)<-3*PI/4)direction[cell*8+0]+=gradientImage.Pixel(f.x+i,f.y+j,0);
                    if (angleImage.Pixel(f.x+i,f.y+j,0)>=-3*PI/4&&angleImage.Pixel(f.x+i,f.y+j,0)<-PI/2)direction[cell*8+1]+=gradientImage.Pixel(f.x+i,f.y+j,0);
                    if (angleImage.Pixel(f.x+i,f.y+j,0)>=-PI/2&&angleImage.Pixel(f.x+i,f.y+j,0)<-PI/4)direction[cell*8+2]+=gradientImage.Pixel(f.x+i,f.y+j,0);
                    if (angleImage.Pixel(f.x+i,f.y+j,0)>=-PI/4&&angleImage.Pixel(f.x+i,f.y+j,0)<0)direction[cell*8+3]+=gradientImage.Pixel(f.x+i,f.y+j,0);
                    if (angleImage.Pixel(f.x+i,f.y+j,0)>=0&&angleImage.Pixel(f.x+i,f.y+j,0)<PI/4)direction[cell*8+4]+=gradientImage.Pixel(f.x+i,f.y+j,0);
                    if (angleImage.Pixel(f.x+i,f.y+j,0)>=PI/4&&angleImage.Pixel(f.x+i,f.y+j,0)<PI/2)direction[cell*8+5]+=gradientImage.Pixel(f.x+i,f.y+j,0);
                    if (angleImage.Pixel(f.x+i,f.y+j,0)>=PI/2&&angleImage.Pixel(f.x+i,f.y+j,0)<3*PI/4)direction[cell*8+6]+=gradientImage.Pixel(f.x+i,f.y+j,0);
                    if (angleImage.Pixel(f.x+i,f.y+j,0)>=3*PI/4&&angleImage.Pixel(f.x+i,f.y+j,0)<PI)direction[cell*8+7]+=gradientImage.Pixel(f.x+i,f.y+j,0);
                }
            }
            
        }
        
        
        for(int cellit=0;cellit<direction.size();cellit++)total_len+=direction[cellit]*direction[cellit];
        total_len=sqrt(total_len);
        for(int cellit=0;cellit<direction.size();cellit++){
            direction[cellit]=direction[cellit]/total_len;
        }
        for(int i=0;i<direction.size();i++)f.data.push_back(direction[i]);
        
        f.data.push_back(Ix.Pixel(f.x, f.y, 0)/gradientImage.Pixel(f.x, f.y, 0));
        f.data.push_back(Iy.Pixel(f.x, f.y, 0)/gradientImage.Pixel(f.x, f.y, 0));
        
        
        
    }
    
    
}


void ComputeMOPSDescriptors(CFloatImage &image, FeatureSet &features)
{
    cout << "Computing mops..." << std::endl;
    // This image represents the window around the feature you need to compute to store as the feature descriptor
    const int windowSize = 8;
    CFloatImage destImage(windowSize, windowSize, 1);
    
    // Gray image
    CFloatImage grayImage = ConvertToGray(image);
    CFloatImage Ix(grayImage.Shape().width,grayImage.Shape().height,1);
    CFloatImage Iy(grayImage.Shape().width,grayImage.Shape().height,1);
    CFloatImage gradientImage(grayImage.Shape().width,grayImage.Shape().height,1);
    CFloatImage angleImage(grayImage.Shape().width,grayImage.Shape().height,1);
    
    
    float Bx[3][3]={{-1,0,1},{-2,0,2},{-1,0,1}};
    float By[3][3]={{-1,-2,-1},{0,0,0},{1,2,1}};
    
    for (int x = 0; x < grayImage.Shape().width; x++) {
        for (int y = 0; y < grayImage.Shape().height; y++) {
            if(x>1&&y>1&&x+1<grayImage.Shape().width-1&&y+1<grayImage.Shape().height-1){
                for (int kx=0;kx<3;kx++){
                    for (int ky=0;ky<3;ky++){
                        Ix.Pixel(x,y,0)=Ix.Pixel(x,y,0)+grayImage.Pixel(x +kx-1,y +ky-1,0)*Bx[3-1-kx][3-1-ky];
                        Iy.Pixel(x,y,0)=Iy.Pixel(x,y,0)+grayImage.Pixel(x+kx-1,y+ky-1,0)*By[3-1-kx][3-1-ky];
                    }
                }
            }
        }
    }
    // Blur Image
    CFloatImage blurImage(image.Shape().width, image.Shape().height, 1);
    
    // Blur 5x5 Gaussian kernel
    int blurWindowSize = 5;
    CFloatImage blurKernel(blurWindowSize, blurWindowSize, 1);
    for(int i = 0; i < blurWindowSize; i++){
        for(int j = 0; j < blurWindowSize; j++){
            blurKernel.Pixel(j,i,0) = gaussian5x5[blurWindowSize * i + j];
        }
    }
    blurKernel.origin[0] = blurWindowSize/2; // x offset to the center of the kernel
    blurKernel.origin[1] = blurWindowSize/2; // y offset to the center of the kernel
    Convolve(grayImage, blurImage, blurKernel);
    
    for (vector<Feature>::iterator i = features.begin(); i != features.end(); i++) {
        Feature &f = *i;
        CTransform3x3 xform;
        
        CTransform3x3 origin;
        origin = origin.Translation( float(-windowSize/2), float(-windowSize/2));
        
        CTransform3x3 rotate;
        rotate = rotate.Rotation(f.angleRadians * 180.0 / PI);
        
        CTransform3x3 trans;
        trans = trans.Translation( f.x, f.y );
        
        CTransform3x3 scale;
        scale[0][0] = 5;
        scale[1][1] = 5;
        
        xform = trans * rotate * scale * origin; //origin * rotate * trans;
        WarpGlobal(blurImage, destImage, xform, eWarpInterpLinear);
        
        f.data.resize(windowSize * windowSize);
        
        // Normalize the patch
        double mean = 0.0;
        for (int row = 0; row < windowSize; row++ ){
            for (int col = 0; col < windowSize; col++){
                mean += destImage.Pixel(col, row, 0);
            }
        }
        mean = mean / (windowSize * windowSize);
        
        
        double variance = 0.0;
        for (int row = 0; row < windowSize; row++ ){
            for (int col = 0; col < windowSize; col++){
                variance += pow((destImage.Pixel(col, row, 0) - mean), 2);
            }
        }
        double stndiv = sqrt( variance / ( (windowSize * windowSize) - 1));
        
        // Fill in the feature descriptor data for a MOPS descriptor
        for (int row = 0; row < windowSize; row++ ){
            for (int col = 0; col < windowSize; col++){
                f.data[row * windowSize + col] = (destImage.Pixel(col, row, 0) - mean)/stndiv;
            }
        }
        
        
        f.data.push_back(0.5*Ix.Pixel(f.x, f.y, 0)/gradientImage.Pixel(f.x, f.y, 0));
        f.data.push_back(0.5*Iy.Pixel(f.x, f.y, 0)/gradientImage.Pixel(f.x, f.y, 0));
        
        
    }
}



//TO DO---------------------------------------------------------------------
//Loop through the image to compute the harris corner values as described in class
// srcImage:  grayscale of original image
// harrisImage:  populate the harris values per pixel in this image
void computeHarrisValues(CFloatImage &srcImage, CFloatImage &harrisImage,CFloatImage &Ix,CFloatImage &Iy)
{
    
    int w = srcImage.Shape().width;//x<w!!!
    int h = srcImage.Shape().height;//y<h!!!
    //CFloatImage Ix(srcImage.Shape().width,srcImage.Shape().height,1);
    //CFloatImage Iy(srcImage.Shape().width,srcImage.Shape().height,1);
    
    float Bx[3][3]={{-1,0,1},{-2,0,2},{-1,0,1}};
    float By[3][3]={{-1,-2,-1},{0,0,0},{1,2,1}};
    
    for (int x = 0; x < w; x++) {
        for (int y = 0; y < h; y++) {
            //printf("\n%f\n",Iy.Pixel(x,y,0));
            if(x>1&&y>1&&x+1<w-1&&y+1<h-1){
                for (int kx=0;kx<3;kx++){
                    for (int ky=0;ky<3;ky++){
                        
                        Ix.Pixel(x,y,0)=Ix.Pixel(x,y,0)+srcImage.Pixel(x +kx-1,y +ky-1,0)*Bx[3-1-kx][3-1-ky];
                        Iy.Pixel(x,y,0)=Iy.Pixel(x,y,0)+srcImage.Pixel(x+kx-1,y+ky-1,0)*By[3-1-kx][3-1-ky];
                        
                    }
                }
            }
        }
    }
    float H[5][5]={{1,4,6,4,1},{4,16,24,16,4},{6,24,36,24,6},{4,16,24,16,4},{1,4,6,4,1}};
    for(int hi=0;hi<5;hi++){
        for(int hj=0;hj<5;hj++){
            H[hi][hj]=(1.0/256.0)*H[hi][hj];
        }
    }
    
    CFloatImage Ix2(srcImage.Shape().width,srcImage.Shape().height,1);
    CFloatImage Iy2(srcImage.Shape().width,srcImage.Shape().height,1);
    CFloatImage Ixy(srcImage.Shape().width,srcImage.Shape().height,1);
    
    for(int x = 0; x < w; x++){
        for (int y = 0; y < h; y++) {
            Ix2.Pixel(x, y, 0)=Ix.Pixel(x,y,0)*Ix.Pixel(x,y,0);
            Iy2.Pixel(x, y, 0)=Iy.Pixel(x,y,0)*Iy.Pixel(x,y,0);
            Ixy.Pixel(x, y, 0)=Ix.Pixel(x,y,0)*Iy.Pixel(x,y,0);
            //angleImage.Pixel(x, y, 0)=atan2(Iy.Pixel(x,y,0), Ix.Pixel(x,y,0));
        }
    }
    
    CFloatImage Ix2_(srcImage.Shape().width,srcImage.Shape().height,1);
    CFloatImage Iy2_(srcImage.Shape().width,srcImage.Shape().height,1);
    CFloatImage Ixy_(srcImage.Shape().width,srcImage.Shape().height,1);
    
    for (int x = 0; x < w; x++) {
        for (int y = 0; y < h; y++) {
            if(x>=2&&y>=2&&x<w-3&&y<h-3){
                for (int kx=0;kx<5;kx++){
                    for (int ky=0;ky<5;ky++){
                        Ix2_.Pixel(x,y,0)+=Ix2.Pixel(x+kx-2,y+ky-2,0)*H[5-1-kx][5-1-ky];
                        Iy2_.Pixel(x,y,0)+=Iy2.Pixel(x+kx-2,y+ky-2,0)*H[5-1-kx][5-1-ky];
                        Ixy_.Pixel(x,y,0)+=Ixy.Pixel(x+kx-2,y+ky-2,0)*H[5-1-kx][5-1-ky];
                    }
                }
            }
        }
    }
    
    for(int x = 0; x < w; x++){
        for (int y = 0; y < h; y++) {
            harrisImage.Pixel(x, y, 0)=(Ix2_.Pixel(x,y,0)*Iy2_.Pixel(x,y,0)-Ixy_.Pixel(x,y,0)*Ixy_.Pixel(x,y,0))/(1e-5+Ix2_.Pixel(x,y,0)+Iy2_.Pixel(x,y,0));
            //printf("\n harrisImage.Pixel(x, y, 0)=%f \n",harrisImage.Pixel(x, y, 0));
        }
    }
    
    
}



// TO DO---------------------------------------------------------------------
// Loop through the harrisImage to threshold and compute the local maxima in a neighborhood
// srcImage:  image with Harris values
// destImage: Assign 1 to a pixel if it is above a threshold and is the local maximum in 3x3 window, 0 otherwise.
//    You'll need to find a good threshold to use.
int computeLocalMaxima(CFloatImage &srcImage,CByteImage &destImage,int r)
{
    int w = srcImage.Shape().width;
    int h = srcImage.Shape().height;
    float threshold=0.05;
    int count=0;
    for (int y = 5; y < h-5; y++) {
        for (int x = 5; x < w-5; x++) {
            destImage.Pixel(x,y,0)=0;
            if( ( x-r>0||y-r>0||x+r<w-1||y+r<h-1))
                if (srcImage.Pixel(x,y,0)>threshold){
                    bool flag=true;
                    for(int dx=-r;dx<=r;dx++){
                        for(int dy=-r;dy<=r;dy++){
                            if( srcImage.Pixel(x, y, 0)<srcImage.Pixel(x+dx, y+dy, 0)&&(dx!=0&&dy!=0)){
                                flag=false;
                            }
                        }
                    }
                    if (flag){
                        count++;
                        destImage.Pixel(x,y,0)=1;
                    }
                }
        }
    }
    return count;
}

// Perform simple feature matching.  This just uses the SSD
// distance between two feature vectors, and matches a feature in the
// first image with the closest feature in the second image.  It can
// match multiple features in the first image to the same feature in
// the second image.
void ssdMatchFeatures(const FeatureSet &f1, const FeatureSet &f2, vector<FeatureMatch> &matches, double &totalScore) {
    printf("computing ssdMatchFeatures...");
    unsigned long m = f1.size();
    unsigned long n = f2.size();
    
    //matches.resize(m);
    FeatureMatch fm;
    totalScore = 0;
    
    double d;
    double dBest;
    int idBest;
    
    for (int i=0; i<m; i++) {
        dBest = 1e100;
        idBest = 0;
        
        for (int j=0; j<n; j++) {
            d = distanceSSD(f1[i].data, f2[j].data);
            
            if (d < dBest) {
                dBest = d;
                idBest = f2[j].id;
            }
        }
        fm.id1=f1[i].id;
        fm.id2=idBest;
        fm.score=dBest;
        //        matches[i].id1 = f1[i].id;
        //        matches[i].id2 = idBest;
        //        matches[i].score = dBest;
        matches.push_back(fm);
        totalScore += fm.score;
        
    }
}

// TODO: Write this function to perform ratio feature matching.
// This just uses the ratio of the SSD distance of the two best matches as the score
// and matches a feature in the first image with the closest feature in the second image.
// It can match multiple features in the first image to the same feature in
// the second image.  (See class notes for more information, and the sshMatchFeatures function above as a reference)
void ratioMatchFeatures(const FeatureSet &f1, const FeatureSet &f2, vector<FeatureMatch> &matches, double &totalScore)
{
    printf("computing ratioMatchFeatures...");
    unsigned long m = f1.size();
    unsigned long n = f2.size();
    
    matches.resize(m);
    totalScore = 0;
    
    //    float d;
    float dBest;
    float dBest2;
    int idBest;
    int idBest2;
    int indexBest=0;
    int indexBest2=0;
    vector<float> matchcost(n);
    
    for (int i=0; i<m; i++) {
        dBest = 1e100;
        idBest = 0;
        dBest2 = 1e100;
        idBest2 = 0;
        //find 1st smallest
        for (int j=0; j<n; j++) {
            matchcost[j] = distanceSSD(f1[i].data, f2[j].data);
            //printf("\n matchcost[j] is %f\n",matchcost[j]);
            if (matchcost[j] < dBest) {
                dBest = matchcost[j];
                idBest = f2[j].id;
                indexBest=j;
            }
        }
        //printf("\n dBest is %f\n",dBest);
        //find 2nd smallest
        for (int j=0; j<n; j++) {
            if (matchcost[j] < dBest2&&j!=indexBest) {
                dBest2 = matchcost[j];
                idBest2 = f2[j].id;
                indexBest2=j;
            }
        }
        //printf("\n dBest2 is %f\n",dBest2);
        matches[i].id1 = f1[i].id;
        matches[i].id2 = idBest;
        matches[i].score = dBest/(dBest2+1e-5);
        totalScore += matches[i].score;
    }
    
}


// Convert Fl_Image to CFloatImage.
bool convertImage(const Fl_Image *image, CFloatImage &convertedImage) {
    if (image == NULL) {
        return false;
    }
    
    // Let's not handle indexed color images.
    if (image->count() != 1) {
        return false;
    }
    
    int w = image->w();
    int h = image->h();
    int d = image->d();
    
    // Get the image data.
    const char *const *data = image->data();
    
    int index = 0;
    
    for (int y=0; y<h; y++) {
        for (int x=0; x<w; x++) {
            if (d < 3) {
                // If there are fewer than 3 channels, just use the
                // first one for all colors.
                convertedImage.Pixel(x,y,0) = ((uchar) data[0][index]) / 255.0f;
                convertedImage.Pixel(x,y,1) = ((uchar) data[0][index]) / 255.0f;
                convertedImage.Pixel(x,y,2) = ((uchar) data[0][index]) / 255.0f;
            }
            else {
                // Otherwise, use the first 3.
                convertedImage.Pixel(x,y,0) = ((uchar) data[0][index]) / 255.0f;
                convertedImage.Pixel(x,y,1) = ((uchar) data[0][index+1]) / 255.0f;
                convertedImage.Pixel(x,y,2) = ((uchar) data[0][index+2]) / 255.0f;
            }
            
            index += d;
        }
    }
    
    return true;
}

// Convert CFloatImage to CByteImage.
void convertToByteImage(CFloatImage &floatImage, CByteImage &byteImage) {
    CShape sh = floatImage.Shape();
    
    assert(floatImage.Shape().nBands == byteImage.Shape().nBands);
    for (int y=0; y<sh.height; y++) {
        for (int x=0; x<sh.width; x++) {
            for (int c=0; c<sh.nBands; c++) {
                float value = floor(255*floatImage.Pixel(x,y,c) + 0.5f);
                
                if (value < byteImage.MinVal()) {
                    value = byteImage.MinVal();
                }
                else if (value > byteImage.MaxVal()) {
                    value = byteImage.MaxVal();
                }
                
                // We have to flip the image and reverse the color
                // channels to get it to come out right.  How silly!
                byteImage.Pixel(x,sh.height-y-1,sh.nBands-c-1) = (uchar) value;
            }
        }
    }
}

// Compute SSD distance between two vectors.
double distanceSSD(const vector<double> &v1, const vector<double> &v2) {
    unsigned long m = v1.size();
    unsigned long n = v2.size();
    
    if (m != n) {
        // Here's a big number.
        return 1e100;
    }
    
    double dist = 0;
    
    for (int i=0; i<m; i++) {
        dist += pow(v1[i]-v2[i], 2);
    }
    
    
    return sqrt(dist);
}

// Transform point by homography.
void applyHomography(double x, double y, double &xNew, double &yNew, double h[9]) {
    double d = h[6]*x + h[7]*y + h[8];
    
    xNew = (h[0]*x + h[1]*y + h[2]) / d;
    yNew = (h[3]*x + h[4]*y + h[5]) / d;
}

// Compute AUC given a ROC curve
double computeAUC(vector<ROCPoint> &results)
{
    double auc=0;
    double xdiff,ydiff;
    for (int i = 1; i < (int) results.size(); i++)
    {
        //fprintf(stream,"%lf\t%lf\t%lf\n",thresholdList[i],results[i].falseRate,results[i].trueRate);
        xdiff=(results[i].falseRate-results[i-1].falseRate);
        ydiff=(results[i].trueRate-results[i-1].trueRate);
        auc=auc+xdiff*results[i-1].trueRate+xdiff*ydiff/2;
        
    }
    return auc;
}