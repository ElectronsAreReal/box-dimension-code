
//1d- und 2d funktionen für praktische benutzbarkeit. zeit- und frequenzbereiche können so als fließkommazahlen eingegeben werden. 
float audioAna::getBoxDimensionSpectrogram(float tmin, float tmax, float fmin, float fmax, bool normalized){
    //returns dimension to user.
    //std::cout << tmin << " " << tmax << " " << fmin << " " << fmax << " " << sR << endl;
    unsigned startPosY=unsigned(tmin*float(sR)/float(step));
    unsigned sizeY=unsigned((tmax-tmin)*float(sR)/float(step));
    unsigned startPosX=unsigned(fmin*float(specF)/sR);
    unsigned stopPosX=unsigned(fmax*float(specF)/sR);
    unsigned sizeX=stopPosX+1-startPosX;
    //std::cout << startPosY << " " << sizeY << " " << startPosX << " " << sizeX << endl;
    if (normalized){
        return getDimension2D(specn,specF,specT,startPosX, startPosY, sizeX, sizeY, float(min(sizeX,sizeY))); 
    }
    else {
        return getDimension2D(spec,specF,specT,startPosX, startPosY, sizeX, sizeY, float(min(sizeX,sizeY)));
    }
}

float audioAna::getDimension2D(float *data2d, unsigned int dataWidth, unsigned int dataHeight, unsigned int startPosX, unsigned int startPosY, unsigned int sizeX, unsigned int sizeY, float strtch){
    unsigned scaleBegin=16;
    unsigned mindim=min(sizeX, sizeY);
    unsigned scalesteps=unsigned(log2(mindim))-unsigned(log2(scaleBegin));
    unsigned stretch=strtch*40.f;
    float boxsizeBegin=float(mindim)/float(scaleBegin);
    float dimension=-1.f;
    if (scalesteps<4){
        cout<<"getDimension2D: size of section to analyze is to small" << endl;
    }
    else{
        float *scales =new float[scalesteps];
        float *boxcounts =new float[scalesteps];
        float *scalesLog =new float[scalesteps];
        float *boxcountsLog =new float[scalesteps];
        boxCount2D(data2d, boxcounts, scales, scalesteps, dataWidth, dataHeight, startPosX, startPosY, sizeX, sizeY, stretch, boxsizeBegin, mindim);
        logarr(scales, scalesLog, scalesteps);
        logarr(boxcounts, boxcountsLog, scalesteps);
        dimension=linReg(scalesLog, boxcountsLog, scalesteps);
    }
    return dimension;
}
//eingentliche boxcount-funktion. 
void audioAna::boxCount2D(float *data, float *bcnt, float *scf, unsigned int scs, unsigned int datW, unsigned int datH, unsigned int startX, unsigned int startY, unsigned int sizeX, unsigned int sizeY, unsigned int stretch, float boxsizeBegin, unsigned int mindim){
    //all dimensions must be power of two
    float upscale=float(stretch);
    float boxs=boxsizeBegin;
    int scsc=0;
    int boxcount;
    float bxu, bxo, byu, byo, bzu, bzo, xb, yb, x, y, z;
    bool allunderbox, alloverbox, openbordercross, top;
    while (boxs>=2){ //loop for each scaling factor
        x=float(startX);
        y=float(startY);
        z=0.;
        boxcount=0;
        while (x<=float(sizeX+startX)-boxs){
            bxu=x;
            bxo=x+boxs;
            y=0.;
            while (y<=float(sizeY+startY)-boxs){
                byu=y;
                byo=y+boxs;
                top=false; //stop lifting the box if all points are under the box
                z=0.;
                while (top==false){
                    if (z>stretch){
                        top=true;
                       // std::cout << " error: maximum height of bosstack reached->check for invalid spectrum data!" << endl;
                    }
                    bzu=z;
                    bzo=z+boxs;
                    allunderbox=true;
                    alloverbox=true;
                    openbordercross=false;
                    xb=bxu;
                    while (xb<bxo && alloverbox==true){//check if all points are over the box
                        yb=byu;
                        while (yb<byo && alloverbox==true){
                            if (data[datW*unsigned(yb)+unsigned(xb)]*upscale<bzo){
                                alloverbox=false;
                            }
                            yb++;
                        }
                        xb++;
                    }

                    if (alloverbox==true){ //if all points are over the box, there cant be one under the box
                        allunderbox=false;
                    }
                    else {
                        xb=bxu;
                        while (xb<bxo && allunderbox==true){//check if all points are under the box
                            yb=byu;
                             while (yb<byo && allunderbox==true){
                                if (data[datW*unsigned(yb)+unsigned(xb)]*upscale>=bzu){
                                    allunderbox=false;
                                }
                                yb++;
                            }
                            xb++;
                        }
                    }

                    if (allunderbox && byo<sizeY+startY && bxo<sizeX+startX){ //every box has three closed sides, where the ponts on the border count, and three open sides where the points do not count. if all points are under or over the box, there is the possibility left that a point at an open border (x or y) is under/over the box. if so, the box needs to be counted, because if one point under the box is next to a point over the box and on the open edge of the box, the plane intersects the box.
                        //if all points are under the box, but at least one point at the open x or y border is higher than the bottom of the box, the plane intersects the box.
                        xb=bxu;
                        while (xb<=bxo && openbordercross==false){
                            if (data[datW*unsigned(byo)+unsigned(xb)]*upscale>=bzu){ //check if one point at the border is higher than the bottom of the box
                                openbordercross=true;
                            }
                            xb++;
                        }
                        yb=byu;
                        while (yb<=byo && openbordercross==false){
                            if (data[datW*unsigned(yb)+unsigned(bxo)]*upscale>=bzu){ //check if one point at the border is higher than the bottom of the box
                                openbordercross=true;
                            }
                            yb++;
                        }
                    }
                    if (alloverbox && byo<sizeY+startY && bxo<sizeX+startX){
                        //if all points are over the box, but at least one point at the open x or y border is lower than the top of the box, the plane intersects the box.
                        xb=bxu;
                        while (xb<=bxo && openbordercross==false){ //check open y-border
                            if (data[datW*unsigned(byo)+unsigned(xb)]*upscale<bzo){ //check if one point at the border is lower than the top of the box
                                openbordercross=true;
                            }
                            xb++;
                        }
                        yb=byu;
                        while (yb<=byo && openbordercross==false){//check open x-border
                            if (data[datW*unsigned(yb)+unsigned(bxo)]*upscale<bzo){ //check if one point at the border is lower than the top of the box
                                openbordercross=true;
                            }
                            yb++;
                        }
                    }
                    if (!(allunderbox || alloverbox) || openbordercross){ //if not all points are under or over the box the box intersects the plane. Or one point at the open x or y border is: lower then the top if all are over the box. or higher than the bottom if all are under the box. Then the plane intersects the box to, even all in the range are under or all are over the box.
                        boxcount++;
                    }

                    if (allunderbox && openbordercross==false){ //if all points are under the box, and the openborderchedk ist negative, the box will not intersect the plane at any position higher than the current one.
                        top=true;
                    }
                    z=z+boxs;
                }
                y=y+boxs;
            }
            x=x+boxs;
        }
        bcnt[scsc]=float(boxcount);
        scf[scsc]=float(mindim)/float(boxs);
        boxs=boxs/2.f;
        //boxs=boxs-2.f;
        scsc++;
    }
}

