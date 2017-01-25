
/////////////////////////////////////////////////////////////////////////////////////////////////
//	Project 4: Eigenfaces                                                                      //
//  CSE 455 Winter 2003                                                                        //
//	Copyright (c) 2003 University of Washington Department of Computer Science and Engineering //
//                                                                                             //
//  File: eigfaces.cpp                                                                         //
//	Author: David Laurence Dewey                                                               //
//	Contact: ddewey@cs.washington.edu                                                          //
//           http://www.cs.washington.edu/homes/ddewey/                                        //
//                                                                                             //
/////////////////////////////////////////////////////////////////////////////////////////////////



#include "stdafx.h"
//#include "FacePosition.h"

EigFaces::EigFaces()
:
Faces()
{
    //empty
}

EigFaces::EigFaces(int count, int width, int height)
:
Faces(count, width, height)
{
    //empty
}

void EigFaces::projectFace(const Face& face, Vector& coefficients) const
{
    if (face.getWidth()!=width || face.getHeight()!=height) {
        throw Error("Project: Face to project has different dimensions");
    }
    
    coefficients.resize(getSize());
    // ----------- TODO #2: compute the coefficients for the face and store in coefficients.
    for(int i=0;i<getSize();i++){
        Face e2;
        face.sub( (*this).getAverage(),e2);
        coefficients[i]= (*this)[i].dot(e2);
    }
}

void EigFaces::constructFace(const Vector& coefficients, Face& result) const
{
    // ----------- TODO #3: construct a face given the coefficients
    result=(*this).getAverage();
    
    for(int pid=0; pid < vector_size; ++pid){
        for(int fid = 0; fid < getSize(); ++fid){
            result[pid] += coefficients[fid] * (*this)[fid][pid];
        }
    }
}

bool EigFaces::isFace(const Face& face, double max_reconstructed_mse,double min_coefficients_mse,  double& mse) const
{
    // ----------- TODO #4: Determine if an image is a face and return true if it is. Return the actual
    // MSE you calculated for the determination in mse
    // Be sure to test this method out with some face images and some non face images
    // to verify it is working correctly.
    Vector coefficients;
    mse = 0;
    this->projectFace(face,coefficients);
    Face result;
    this->constructFace(coefficients, result);
    for(int i=0;i<face.getSize();i++){
        mse += (result[i]-face[i]) *  (result[i]-face[i]);
    }
    mse /= (double)face.getSize();
    return (mse<max_reconstructed_mse)&&(mse>min_coefficients_mse);  // placeholder, replace
}

bool EigFaces::verifyFace(const Face& face, const Vector& user_coefficients, double max_coefficients_mse,double& mse) const
{
    // ----------- TODO #5 : Determine if face is the same user give the user's coefficients.
    // return the MSE you calculated for the determination in mse.
    Vector coefficients;
    this->projectFace(face,coefficients);
    for(int i=0;i<coefficients.getSize();i++){
        mse+=pow((user_coefficients[i]-coefficients[i]),2);
    }
    mse/=coefficients.getSize();
    return (mse<max_coefficients_mse);  // placeholder, replace
}

void EigFaces::recognizeFace(const Face& face, Users& users) const
{
    // ----------- TODO #6: Sort the users by closeness of match to the face
    
    Vector coefficients;
    Vector mse;
    this->projectFace(face,coefficients);
    for(int j=0;j<users.getSize();j++){
        users[j].sub(coefficients, mse);
        users[j].setMse(pow(mse.mag(),2)/coefficients.getSize());
    }
    users.sort();
}



void EigFaces::findFace(const Image& img, double min_scale, double max_scale, double step, int n, bool crop, Image& result) const
{
    // ----------- TODO #7: Find the faces in Image. Search image scales from min_scale to max_scale inclusive,
    // stepping by step in between. Find the best n faces that do not overlap each other. If crop is true,
    // n is one and you should return the cropped original img in result. The result must be identical
    // to the original besides being cropped. It cannot be scaled and it must be full color. If crop is
    // false, draw green boxes (use r=100, g=255, b=100) around the n faces found. The result must be
    // identical to the original image except for the addition of the boxes.
    Image grayimage;
    img.grayscale(grayimage);
    
    auto isOverlap = [&](const FacePosition& f1, const FacePosition& f2){
        const int fw1 = getWidth() / f1.scale;
        const int fh1 = getHeight() / f1.scale;
        const int fw2 = getWidth() / f2.scale;
        const int fh2 = getHeight() / f2.scale;
        
        if(f1.x / f1.scale > f2.x / f2.scale + fw2 || f1.x / f1.scale + fw1 < f2.x / f2.scale ||
           f1.y / f1.scale > f2.y / f2.scale + fh2 || f1.y / f1.scale + fh1 < f2.y / f2.scale){
            return false;
        }
        return true;
    };
        
    printf("Number of eigen face: %d\n", getSize());
    
    std::vector<FacePosition> list_faceposition;
    list_faceposition.reserve(n);
    for(double scale = min_scale; scale <= max_scale; scale += step){
        printf("Searching at scale %.3f\n", scale);
        Face sample_face(img.getWidth()*scale,img.getHeight()*scale);
        
        grayimage.resample(sample_face);
        printf("Resampled size: (%d,%d)\n", sample_face.getWidth(), sample_face.getHeight());
        for(int dx=0;dx<sample_face.getWidth()-getWidth(); ++dx){
            for(int dy=0;dy<sample_face.getHeight()-getHeight(); ++dy){
                Face box(getWidth(),getHeight());
                Face gradientbox(getWidth(),getHeight());
                box.subimage(dx, dx+getWidth()-1, dy, dy+getHeight()-1, sample_face, 1);
                FacePosition boxp;
                
                
                if(isFace(box, 2000, 200, boxp.error)){//use 200-2000 to store smaller number of candidate boxes
                    boxp.x=dx;
                    boxp.y=dy;
                    boxp.scale=scale;
                    Face distance_2_mean;
                    box.sub(this->getAverage(), distance_2_mean);
                    boxp.error=boxp.error*distance_2_mean.mag()/(box.var());
                    list_faceposition.push_back(boxp);
                }
            }
        }
    }
    
    std::sort(list_faceposition.begin(), list_faceposition.end(),
              [](const FacePosition& b1, const FacePosition& b2){return b1.error < b2.error;});
    
    if(list_faceposition.empty()){
        printf("list is empty\n");
        exit(-1);
    }
    
    vector<bool> is_valid(list_faceposition.size(), true);
    
    result = img;
    int index = 0;
    for(int i=0; i<list_faceposition.size(); ++i){
        const FacePosition& m = list_faceposition[i];
        if(index == n){
            break;
        }
        
        for(int j=0; j<i; ++j){
            if(is_valid[j]){
                if(isOverlap(list_faceposition[i], list_faceposition[j])){
                    is_valid[i] = false;
                }
            }
        }
        
        if(!is_valid[i]){
            continue;
        }
        
        const double& scale = m.scale;
        
        int xl = m.x / scale;
        int yl = m.y / scale;
        int xh = floor((m.x + getWidth()) / scale);
        int yh = floor((m.y + getHeight()) / scale);
        Face result_check_r(xh-xl+1,yh-yl+1);
        Face result_check_g(xh-xl+1,yh-yl+1);
        Face result_check_b(xh-xl+1,yh-yl+1);
        double divcolor1 = 0.0;
        double divcolor2 = 0.0;
        for(int i1=xl;i1<xh;i1++){
            for(int i2=yl;i2<yh;i2++){
                result_check_r.pixel(i1-xl, i2-yl, 0)=img.pixel(i1,i2,0);
                result_check_g.pixel(i1-xl, i2-yl, 0)=img.pixel(i1,i2,1);
                result_check_b.pixel(i1-xl, i2-yl, 0)=img.pixel(i1,i2,2);
            }
        }
        divcolor1=result_check_r.mag()-result_check_g.mag();
        divcolor2=result_check_r.mag()-result_check_b.mag();
        //we think a face is closer to red than blue and green
        if(divcolor1<0)continue;
        if(divcolor2<0)continue;
        // only for 'crop my own photo'
        if(m.x/scale<(xh-xl)/10)continue;
        
        std::cout<<m.x/scale<<" "<<m.y/scale<<" "<<m.scale<<" "<<m.error<<" divcolor is : "<<divcolor1<<" "<<divcolor2<<std::endl;
        if(!crop){
        result.line(xl, yl, xl,yh, 100, 255, 100);
        result.line(xl,yl, xh, yl, 100, 255, 100);
        result.line(xh,yl, xh, yh, 100, 255, 100);
        result.line(xl,yh, xh, yh, 100, 255, 100);
        }
        else{
        img.crop(xl,yl,xh,yh,result);
        }
        index++;
    }
    
}

void EigFaces::morphFaces(const Face& face1, const Face& face2, double distance, Face& result) const
{
    // TODO (extra credit): MORPH along *distance* fraction of the vector from face1 to face2 by
    // interpolating between the coefficients for the two faces and reconstructing the result.
    // For example, distance 0.0 will approximate the first, while distance 1.0 will approximate the second.
    // Negative distances are ok two.
    
}

const Face& EigFaces::getAverage() const
{
    return average_face;
}

void EigFaces::setAverage(const Face& average)
{
    average_face=average;
}



