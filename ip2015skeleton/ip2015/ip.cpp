#include "ip.h"
#include "main.h"
#include <algorithm>
#include <stdlib.h>
#include <cmath>
#include <time.h>
#include <stack>



/*
 * convolve with a box filter
 */
Image* ip_blur_box (Image* src, int size)
{
    cerr << "This function is not implemented." << endl;
    return NULL;
}

/*
 * convolve with a gaussian filter
 */
Image* ip_blur_gaussian (Image* src, int size, double sigma)
{
    cerr << "This function is not implemented." << endl;
    return NULL;
}


/*
 * convolve with a triangle filter
 */
Image* ip_blur_triangle (Image* src, int size)
{
    cerr << "This function is not implemented." << endl;
    return NULL;
}


/*
 * interpolate with a black image
 */
Image* ip_brighten (Image* src, double alpha)
{
    cerr << "This function is not implemented." << endl;
    return NULL;
}


/*
 * shift colors
 */
Image* ip_color_shift(Image* src)
{
    cerr << "This function is not implemented." << endl;
    return NULL;
}


/*
 * use a mask image for a per-pixel alpha value to perform
 * interpolation with a second image
 */
Image* ip_composite (Image* src1, Image* src2,
                     Image* mask)
{
    cerr << "This function is not implemented." << endl;
    return NULL;
}


/*
 * interpolate with the average intensity of the src image
 */
Image* ip_contrast (Image* src, double alpha)
{
    cerr << "This function is not implemented." << endl;
    return NULL;
}



/*
 * convolve an image with a kernel
 */
Image* ip_convolve (Image* src, int size, double* kernel )
{
    cerr << "This function is not implemented." << endl;
    return NULL;
}



/*
 *  create cropped version of image
 */
Image* ip_crop (Image* src, int x0, int y0, int x1, int y1)
{
    cerr << "This function is not implemented." << endl;
    return NULL;
}
/*
 * convolve with an edge detection kernel
 */
Image* ip_edge_detect (Image* src)
{
    cerr << "This function is not implemented." << endl;
    return NULL;
}


/*
 * extract channel of input image
 */
Image* ip_extract (Image* src, int channel)
{
    cerr << "This function is not implemented." << endl;
    return NULL;
}


/*
 * create your own fun warp
 */
//Image* ip_fun_warp (Image* src)
//{
//    cerr << "This function is not implemented." << endl;
//    //  ask user for input parameters here including resampling method and,
//    //  if gaussian resampling is used, its filtersize and sigma
//    //  if you implement more than one warp, you should ask the
//    //  user to chose the one to perform here too!
//    
//    return NULL;
//}
/*
 * create a new image with values equal to the psychosomatic intensities
 * of the source image
 */
Image* ip_grey (Image* src)
{
    cerr << "This function is not implemented." << endl;
    return NULL;
}


/*
 *  shift image by dx and dy (modulo width & height)
 */

Image* ip_image_shift (Image* src, double dx, double dy)
{
    cerr << "This function is not implemented." << endl;
    return NULL;
}
/*
 * interpolate an image with another image
 */
Image* ip_interpolate (Image* src1, Image* src2, double alpha)
{
    cerr << "This function is not implemented." << endl;
    return NULL;
}
/*
 * invert input image
 */
Image* ip_invert (Image* src)
{
    cerr << "This function is not implemented." << endl;
    return NULL;
}
bool comparePixel(Pixel &pix1, Pixel &pix2){
    bool similar = true;
    for (int i = 0; i < 3; ++i) similar &= (abs(pix1.getColor(i) - pix2.getColor(i)) < .05);
    return similar;
}

void checkNeighbors(int x, int y, Image &src, vector<bool>* neighbors){
    int srcWidth = src.getWidth();
    int srcHeight = src.getHeight();
    
    Pixel pix = src.getPixel(x, y);
    Pixel pixN;
    
    if(x > 0){
        if (y > 0) {
            pixN = src.getPixel(x-1,y-1);
            neighbors->at(0) = comparePixel(pix,pixN);
        }
        pixN = src.getPixel(x-1,y);
        neighbors->at(3) = comparePixel(pix, pixN);
        if (y < srcHeight - 1) {
            pixN = src.getPixel(x-1,y+1);
            neighbors->at(5) = comparePixel(pix, pixN);
        }
    }
    if (x < srcWidth - 1) {
        if (y > 0) neighbors->at(2) = comparePixel(pix, src.getPixel(x+1, y-1, pixN));
        neighbors->at(4) = comparePixel(pix, src.getPixel(x+1, y, pixN));
        if (y < srcHeight - 1) neighbors->at(7) = comparePixel(pix, src.getPixel(x+1, y+1, pixN));
    }
    if(y > 0) neighbors->at(1) = comparePixel(pix, src.getPixel(x, y-1, pixN));
    if(y < srcHeight -1) neighbors->at(6) = comparePixel(pix, src.getPixel(x, y+1, pixN));
     
    
}

void drawEdge(int startX, int startY, int endX, int endY, Image &src){
    Pixel black(0,0,0);
    if(endX > startX){
        if(endY > startY){
            for(int i = 0; i < endX - startX; i++){
                src.setPixel(startX+i, startY+i, black);
            }
        }else{
            for(int i = 0; i<endX - startX; i++){
                src.setPixel(startX+i, startY-i, black);
            }
        }
    }else{
        if(endY > startY){
            for(int i = 0; i < endX - startX; i++){
                src.setPixel(startX-i, startY+i, black);
            }
        }else{
            for(int i = 0; i<endX - startX; i++){
                src.setPixel(startX-i, startY-i, black);
            }
        }
    }
}

void drawLine(int startX, int startY, int endX, int endY, int blockSize, Image &src){
    bool xChanges = startX != endX;
    bool yChanges = startY != endY;
    Pixel red(1,0,0);
    if(xChanges){
        if(yChanges){
            if(endX > startX){
                if(endY > startY){
                    for(int i = 0; i < endX - startX; i++){
                        src.setPixel(startX+i, startY+i, red);
                    }
                }else{
                    for(int i = 0; i<endX - startX; i++){
                        src.setPixel(startX+i, startY-i, red);
                    }
                }
            }else{
                if(endY > startY){
                    for(int i = 0; i < endX - startX; i++){
                        src.setPixel(startX-i, startY+i, red);
                    }
                }else{
                    for(int i = 0; i<endX - startX; i++){
                        src.setPixel(startX-i, startY-i, red);
                    }
                }
            }
        }else{
            int i = min(startX, endX);
            int bigI = max(startX, endX);
            int j = startY;
            while(i <= bigI){
                src.setPixel(i, j, red);
                i++;
            }
        }
    }else if(yChanges){
        int j = min(startY, endY);
        int bigJ = max(startY, endY);
        int i = startX;
        while(j <= bigJ){
            src.setPixel(i, j, red);
            j++;
        }
    }
}



void drawEdge(int startX, int startY, int endX, int endY, int blockSize, Image &src){
    int startCenterX = startX*blockSize + blockSize/2;
    int startCenterY = startY*blockSize + blockSize/2;
    int endCenterY = endY*blockSize + blockSize/2;
    int endCenterX = endX*blockSize + blockSize/2;
    drawLine(startCenterX,startCenterY,endCenterX,endCenterY,blockSize,src);
}

void getNeighbor(int startX, int startY, int* neighborX, int* neighborY, int index){
    switch(index){
        case 0: {
            *neighborX = startX-1;
            *neighborY = startY-1;
            break;
        }
        case 1: {
            *neighborX = startX;
            *neighborY = startY-1;
            break;
        }
        case 2: {
            *neighborX = startX+1;
            *neighborY = startY-1;
            break;
        }
        case 3: {
            *neighborX = startX-1;
            *neighborY = startY;
            break;
        }
        case 4: {
            *neighborX = startX+1;
            *neighborY = startY;
            break;
        }
        case 5: {
            *neighborX = startX-1;
            *neighborY = startY+1;
            break;
        }
        case 6: {
            *neighborX = startX;
            *neighborY = startY+1;
            break;
        }
        case 7: {
            *neighborX = startX+1;
            *neighborY = startY+1;
            break;
        }
    }
}

int getNeighborIndex(int startIndex, int blockSizeX, int blockSizeY, int neighborIndex){
    int resultIndex = -1;
    switch(neighborIndex){
        case 0: {
            resultIndex = startIndex - blockSizeY - 1;
            break;
        }
        case 1: {
            resultIndex = startIndex - 1;
            break;
        }
        case 2: {
            resultIndex = startIndex + blockSizeY - 1;
            break;
        }
        case 3: {
            resultIndex = startIndex - blockSizeY;
            break;
        }
        case 4: {
            resultIndex = startIndex + blockSizeY;
            break;
        }
        case 5: {
            resultIndex = startIndex - blockSizeY + 1;
            break;
        }
        case 6: {
            resultIndex = startIndex + 1;
            break;
        }
        case 7: {
            resultIndex = startIndex + blockSizeY + 1;
            break;
        }
    }
    return resultIndex;
}

// Pixel at i,j represents left corner of 2x2 block to check
void checkRedundancy(int i,
                     int j,
                     int blockSizeX,
                     int blockSizeY,
                     vector< vector<bool> > * similarityGraph) {
    vector<bool> adjList = similarityGraph->at(i*blockSizeY + j);
    bool upperLeft = adjList[4] && adjList[6] && adjList[7];
    adjList = similarityGraph->at(i*blockSizeY + j + 1);
    bool bottomLeft = adjList[2] && adjList[4];
    adjList = similarityGraph->at((i+1)*blockSizeY + j + 1);
    bool bottomRight = adjList[1];
    if (upperLeft && bottomRight && bottomLeft) {
        similarityGraph->at(i*blockSizeY + j)[7] = false;
        similarityGraph->at((i+1)*blockSizeY + j+1)[0] = false;
        similarityGraph->at(i*blockSizeY + j+1)[2] = false;
        similarityGraph->at((i+1)*blockSizeY + j)[5] = false;
    }
}

/*
 * Borders is: 
 * {left, right, top, bottom}
 */
bool bounds(int i,
            int j,
            int * borders) {
    return (borders[0] <= i && i <= borders[1] && borders[2] <= j && j <= borders[3]);
}

/*
 * Returns the size of the connected component at pixel i, j.
 */
int sizeConnectedComponent(int i,
                           int j,
                           int blockSizeX,
                           int blockSizeY,
                           int * borders,
                           vector< vector<bool> > * similarityGraph) {
    vector<bool> visited(blockSizeX*blockSizeY, false);
    stack<int> nodesToVisit;
    vector<bool> adjlist = similarityGraph->at(i*blockSizeY+j);
    int neighborX = i;
    int neighborY = j;
    int numberOfConnectedComponents = 0;
    for (int index = 0; index < adjlist.size(); ++index) {
        if(adjlist[index] && bounds(i, j, borders)) {
            getNeighbor(i, j, &neighborX, &neighborY , index);
            nodesToVisit.push(neighborX * blockSizeY + neighborY);
            visited[neighborX * blockSizeY + neighborY] = true;
            numberOfConnectedComponents++;
        }
    }

    int neighbor = i * blockSizeY + blockSizeX;
    while (!nodesToVisit.empty()) {
        int visitedNode = nodesToVisit.top();
        adjlist = similarityGraph->at(visitedNode);
        for (int index = 0; index < adjlist.size(); ++index) {
            getNeighbor(visitedNode / blockSizeY, visitedNode % blockSizeY, &neighborX, &neighborY, index);
            neighbor = neighborX * blockSizeY + neighborY;
            if (adjlist[index] && !visited[neighbor] && bounds(neighborX, neighborY, borders)) {
                nodesToVisit.push(neighbor);
                visited[neighbor] = true;
                numberOfConnectedComponents++;
            }
            
        }
        nodesToVisit.pop();
    }
    
    return numberOfConnectedComponents;
}

void sparsePixels(int i,
                  int j,
                  int blockSizeX,
                  int blockSizeY,
                  int * weightMajor,
                  int * weightMinor,
                  vector< vector<bool> > * similarityGraph) {
    int left = max(0, i - 3);
    int right = min(blockSizeX, i + 4);
    int top = max(0, j - 3);
    int bottom = min(blockSizeY, j + 4);
    int borders[] = {left, right, top, bottom};
    int major = sizeConnectedComponent(i + 1, j, blockSizeX, blockSizeY, borders, similarityGraph);
    int minor = sizeConnectedComponent(i, j, blockSizeX, blockSizeY, borders, similarityGraph);
    
    if (major > minor) *weightMinor += major - minor;
    else *weightMajor += minor - major;
}

int getNodeDegree(vector<bool> & neighbors) {
    int total = 0;
    for (int i = 0; i < neighbors.size(); ++i) {
        if (neighbors[i]) ++total;
    }
    return total;
}

/*
 * Check if one of the diagonal nodes are of valence one. If so, add 5 to
 * its weight.
 */
void islands(int i,
             int j,
             int blockSizeX,
             int blockSizeY,
             int * weightMajor,
             int * weightMinor,
             vector< vector<bool> > & similarityGraph) {
    int topLeft = getNodeDegree(similarityGraph[i*blockSizeY+j]);
    int topRight = getNodeDegree(similarityGraph[i*blockSizeY+j+1]);
    int botLeft = getNodeDegree(similarityGraph[(i+1)*blockSizeY+j]);
    int botRight =getNodeDegree(similarityGraph[(i+1)*blockSizeY+j+1]);
    if (topLeft == 1 || botRight == 1) weightMinor += 5;
    //question
    if (topRight == 1 || botLeft == 1) weightMajor += 5;
}

/*
 *goes to the next node given a node index and a visited path. Presuming that the node in question is a valence-2 node
 */
int goToNextNode(int nodeIndex, int blockSizeX, int blockSizeY, int visited, vector< vector<bool> > & similarityGraph){
    vector<bool> neighbors = similarityGraph[nodeIndex];
    for(int i = 0; i < neighbors.size(); i++){
        if(neighbors[i]){
            int neighborIndex = getNeighborIndex(nodeIndex, blockSizeX, blockSizeY, i);
            if(neighborIndex != visited){
                return neighborIndex;
            }
        }
    }
    return -1;
}

void curve(int i,
           int j,
           int blockSizeX,
           int blockSizeY,
           int * weightMajor,
           int * weightMinor,
           vector< vector<bool> > & similarityGraph) {
    int topLeftIndex = i*blockSizeY+j;
    int topRightIndex = i*blockSizeY+j+1;
    int botLeftIndex = (i+1)*blockSizeY+j;
    int botRightIndex = (i+1)*blockSizeY+j+1;
    int topLeft = getNodeDegree(similarityGraph[topLeftIndex]);
    int topRight = getNodeDegree(similarityGraph[topRightIndex]);
    int botLeft = getNodeDegree(similarityGraph[botLeftIndex]);
    int botRight =getNodeDegree(similarityGraph[botRightIndex]);
    int topLeftVisited = botRightIndex;
    int topRightVisited = botLeftIndex;
    int botLeftVisited = topRightIndex;
    int botRightVisited = topLeftIndex;
    int newTopLeftIndex, newTopRightIndex, newBotLeftIndex, newBotRightIndex;
    int majorLength = 1;
    int minorLength = 1;
    while(topLeft == 2){
        minorLength++;
        newTopLeftIndex = goToNextNode(topLeftIndex, blockSizeX, blockSizeY, topLeftVisited, similarityGraph);
        topLeftVisited = topLeftIndex;
        topLeftIndex = newTopLeftIndex;
        topLeft = getNodeDegree(similarityGraph[newTopLeftIndex]);
    }
    while(botRight == 2){
        minorLength++;
        newBotRightIndex = goToNextNode(botRightIndex, blockSizeX, blockSizeY, botRightVisited, similarityGraph);
        botRightVisited = botRightIndex;
        botRightIndex = newBotRightIndex;
        botRight = getNodeDegree(similarityGraph[newBotRightIndex]);
    }
    while(topRight == 2){
        majorLength++;
        newTopRightIndex = goToNextNode(topRightIndex, blockSizeX, blockSizeY, topRightVisited, similarityGraph);
        topRightVisited = topRightIndex;
        topRightIndex = newTopRightIndex;
        topRight = getNodeDegree(similarityGraph[newTopRightIndex]);
    }
    while(botLeft == 2){
        majorLength++;
        newBotLeftIndex = goToNextNode(botLeftIndex, blockSizeX, blockSizeY, botLeftVisited, similarityGraph);
        botLeftVisited = botLeftIndex;
        botLeftIndex = newBotLeftIndex;
        botLeft = getNodeDegree(similarityGraph[newBotLeftIndex]);
    }
    if(majorLength>minorLength){
        weightMajor += 6;
    }else if(minorLength>majorLength){
        weightMinor += 6;
    }
}

/*
 * Pixel at i,j is the left hand corner of a 2x2 block with diagonals.
 */
void chooseDiagonals(int i,
                     int j,
                     int blockSizeX,
                     int blockSizeY,
                     vector< vector<bool> > * similarityGraph) {
    int weightMajor = 0;
    int weightMinor = 0;
    sparsePixels(i, j, blockSizeX, blockSizeY, &weightMajor, &weightMinor, similarityGraph);
    islands(i, j, blockSizeX, blockSizeY, &weightMajor, &weightMinor, *similarityGraph);
    curve(i, j, blockSizeX, blockSizeY, &weightMajor, &weightMinor, *similarityGraph);
    if (weightMinor > weightMajor) {
        similarityGraph->at((i)*blockSizeY+j+1)[2] = false;
        similarityGraph->at((i+1)*blockSizeY+j)[5] = false;
    } else {
        similarityGraph->at((i)*blockSizeY+j)[7] = false;
        similarityGraph->at((i+1)*blockSizeY+j+1)[0] = false;
    }
}

void drawLineFromMidpoint(int i, int j, int pixelSize, int neighborIndex, Image& src, Pixel p){
    int halfRange = pixelSize/4;
    if(neighborIndex == 0 || neighborIndex == 7){
        for(int k = - halfRange; k <= halfRange; k++){
            src.setPixel(i + k, j - k, p);
            src.setPixel(i + k + 1, j - k, p);
            src.setPixel(i + k - 1, j - k, p);
        }
    }else if(neighborIndex == 2 || neighborIndex == 5){
        for(int k = - halfRange ; k <= halfRange; k++){
            src.setPixel(i + k, j + k, p);
            src.setPixel(i + k - 1, j + k, p);
            src.setPixel(i + k + 1, j + k, p);
        }
    }
}

void reshapePixel(int i, int j, int pixelSize, int neighborIndex, Image& src){
    
    Pixel p = src.getPixel(i*pixelSize + pixelSize/2+1, j*pixelSize + pixelSize/2+2);
    int endDiag;
    if(neighborIndex == 0){
        int startDiagX = i*pixelSize;
        int startDiagY = j*pixelSize;
        for(int k = 0; k <= pixelSize/4+1; k++){
            drawLineFromMidpoint(startDiagX+k, startDiagY+k, pixelSize, 0, src, p);
        }
    }else if(neighborIndex == 2){
        int startDiagX = i*pixelSize + pixelSize;
        int startDiagY = j*pixelSize;
        for(int k = 0; k <= pixelSize/4+1; k++){
            drawLineFromMidpoint(startDiagX-k, startDiagY+k, pixelSize, 2, src, p);
        }
    }else if(neighborIndex == 5){
        int startDiagX = i*pixelSize;
        int startDiagY = j*pixelSize + pixelSize;
        for(int k = 0; k <= pixelSize/4+1; k++){
            drawLineFromMidpoint(startDiagX+k, startDiagY-k, pixelSize, 5, src, p);
        }
    }else if(neighborIndex == 7){
        int startDiagX = i*pixelSize + pixelSize;
        int startDiagY = j*pixelSize + pixelSize;
        for(int k = 0; k <= pixelSize/4+1; k++){
            drawLineFromMidpoint(startDiagX-k, startDiagY-k, pixelSize, 7, src, p);
        }
    }
}

void reshapePixels(vector< vector<bool> > &similarityGraph,
                 Image& src,
                 int blockSizeX,
                 int blockSizeY,
                 int pixelSize){
    for(int i = 0; i < blockSizeX; i++){
        for(int j = 0; j < blockSizeY; j++){
            vector<bool> adjlist = similarityGraph[i*blockSizeY+j];
            for (int k = 0; k < adjlist.size(); ++k)
                if (adjlist[k]) reshapePixel(i, j, pixelSize, k, src);
        }
    }
}

void reshapePixelGl(int i,
                    int j,
                    vector<bool> adjList,
                    ImageGraph* src,
                    int blockSizeX,
                    int blockSizeY,
                    int pixelSize){

    ImagePixel* p = &currentImageGraph->at(i*blockSizeY+j);
    //cout << p << endl;

    ImagePixel* T;
    ImagePixel* TL;
    ImagePixel* TR;
    ImagePixel* L;
    ImagePixel* R;
    
    if(adjList[0]){
        TL = &(*src)[(i-1)*blockSizeY + (j-1)];
        T = &(*src)[i*blockSizeY + (j-1)];
        L = &(*src)[(i-1)*blockSizeY + j];
        
        p->vertices[0][0] += pixelSize/4;
        p->vertices[0][1] += pixelSize/4;
        p->vertices[1][0] -= pixelSize/4;
        p->vertices[1][1] -= pixelSize/4;
       // cout << p->vertices[0][0] << "," << p -> vertices[0][1] << endl;
        
        
        T->vertices[2][0] += 0.25*pixelSize;
        T->vertices[2][1] += 0.25*pixelSize;
        T->vertices[3][0] += 0.25*pixelSize;
        T->vertices[3][1] += 0.25*pixelSize;
        
        TL->vertices[4][0] -= 0.25*pixelSize;
        TL->vertices[4][1] -= 0.25*pixelSize;
        TL->vertices[5][0] += 0.25*pixelSize;
        TL->vertices[5][1] += 0.25*pixelSize;
        
        L->vertices[6][0] -= 0.25*pixelSize;
        L->vertices[6][1] -= 0.25*pixelSize;
        L->vertices[7][0] -= 0.25*pixelSize;
        L->vertices[7][1] -= 0.25*pixelSize;
    }
    
    if(adjList[2]){
        T = &(*src)[i*blockSizeY + (j-1)];
        TR = &(*src)[(i+1)*blockSizeY + (j-1)];
        R = &(*src)[(i+1)*blockSizeY + j];
        
        p->vertices[6][0] += 0.25*pixelSize;
        p->vertices[6][1] -= 0.25*pixelSize;
        p->vertices[7][0] -= 0.25*pixelSize;
        p->vertices[7][1] += 0.25*pixelSize;
//        

        T->vertices[4][0] -= 0.25*pixelSize;
        T->vertices[4][1] += 0.25*pixelSize;
        T->vertices[5][0] -= 0.25*pixelSize;
        T->vertices[5][1] += 0.25*pixelSize;
        
        TR->vertices[2][0] -= 0.25*pixelSize;
        TR->vertices[2][1] += 0.25*pixelSize;
        TR->vertices[3][0] += 0.25*pixelSize;
        TR->vertices[3][1] -= 0.25*pixelSize;
        
        R->vertices[0][0] += 0.25*pixelSize;
        R->vertices[0][1] -= 0.25*pixelSize;
        R->vertices[1][0] += 0.25*pixelSize;
        R->vertices[1][1] -= 0.25*pixelSize;
    }
    
//
//    ImagePixel TR = src[(i+1)*blockSizeY + (j-1)];
    
}

void reshapePixelsGl(vector< vector<bool> > &similarityGraph,
                     ImageGraph* src,
                     int blockSizeX,
                     int blockSizeY,
                     int pixelSize){
    for(int i=0; i<blockSizeX; i++){
        for(int j=0; j<blockSizeY; j++){
            reshapePixelGl(i, j, similarityGraph[i*blockSizeY+j], src, blockSizeX, blockSizeY, pixelSize);
        }
    }

}



/*
 * define your own filter
 * you need to request any input parameters here, not in control.cpp
 */

Image* ip_misc(Image* src,
               ImageGraph* imageGraph,
               vector<vector<bool>>* similarityGraph,
               const int blockSizeX,
               const int blockSizeY,
               const int pixelSize)
{
    //cerr << "This function is not implemented." << endl;
//    const int blockSizeY = 16;
//    const int blockSizeX = 40;
//    const int pixelSize = 15;
//
    int srcWidth = src->getWidth();
    int srcHeight = src->getHeight();
    Image* rawGraph = new Image(blockSizeX, blockSizeY, 8);
    currentImageGraph = new ImageGraph(blockSizeX * blockSizeY);
    
    Pixel srcPixel;
    
    int blockWidth = srcWidth / blockSizeX;
    int blockHeight = srcHeight / blockSizeY;
    
    for(int i = 0; i<blockSizeX; i++){
        for(int j=0; j<blockSizeY; j++){
            Pixel outputPixel(0,0,0);
            int blockEndi = blockWidth * i + blockWidth;
            int blockEndj = blockHeight * j + blockHeight;
            blockEndi = min(blockEndi, srcWidth);
            blockEndj = min(blockEndj, srcHeight);
            outputPixel =
              src->getPixel(i * blockWidth + blockWidth / 2, j * blockHeight + blockHeight / 2);
            
            vector<vector<GLfloat>> vertices(8, vector<GLfloat>(2));
            
            // Polygon ordering counter clockwise
            
            for (int k = 0; k < 4; ++k) {
                vertices[2*k][0] = static_cast<float>((i + (k > 1 ? 1: 0)) * pixelSize);
                vertices[2*k][1] = static_cast<float>(16 * pixelSize - (j + (k==1 || k ==2 ? 1:0)) * pixelSize);
                vertices[2*k+1][0] = static_cast<float>((i + (k > 1 ? 1: 0)) * pixelSize);
                vertices[2*k+1][1] = static_cast<float>(16 * pixelSize - (j + (k==1 || k ==2 ? 1:0)) * pixelSize);
            }
            
            currentImageGraph->at(i*blockSizeY+j) = ImagePixel(outputPixel, vertices);
            rawGraph->setPixel(i, j, outputPixel);
//            outputPixel = Pixel(0,0,0);
//            rawGraph->setPixel(i, j, outputPixel);
            
        }
    }
    
    similarityGraph = new vector<vector<bool>>(blockSizeX * blockSizeY, vector<bool>(8,false));
    
    for (int i = 0; i < blockSizeX; ++i) {
        for (int j = 0; j < blockSizeY; ++j) {
            checkNeighbors(i, j, *rawGraph, &similarityGraph->at(i*blockSizeY + j));
        }
    }
    
//    
//    for (int i = 0; i < blockSizeX; ++i) {
//        for (int j = 0; j < blockSizeY; ++j) {
//            for (int k = 0; k < 8; ++k) cout << similarityGraph->at(i*blockSizeY + j)
//                [k] << ",";
//            cout << endl;
//        }
//    }
//    
//    
    
    
    Image* rawGraphTest = new Image(blockSizeX*pixelSize, blockSizeY*pixelSize, 3);
    
    Pixel rawPixel(0,0,0);
    
    for(int i=0; i<blockSizeX; ++i){
        for(int j=0; j<blockSizeY; ++j){
            rawPixel = rawGraph->getPixel(i, j);
            for(int x = i*pixelSize; x < i*pixelSize+pixelSize; ++x){
                for(int y = j * pixelSize; y < j*pixelSize + pixelSize; ++y){
                    rawGraphTest->setPixel(x, y, rawPixel);
                }
            }
        }
    }

    for (int i = 0; i < blockSizeX - 1; ++i) {
        for (int j = 0; j < blockSizeY - 1; ++j) {
            // If there are cross edges
            if (similarityGraph->at(i*blockSizeY + j)[7] && similarityGraph->at(i*blockSizeY + j + 1)[2]) {
                checkRedundancy(i,j,blockSizeX, blockSizeY, similarityGraph);
            }
            
        }
    }
    for (int i = 0; i < blockSizeX - 1; ++i) {
        for (int j = 0; j < blockSizeY - 1; ++j) {
            if (similarityGraph->at(i*blockSizeY + j)[7] && similarityGraph->at(i*blockSizeY + j + 1)[2]) {
                chooseDiagonals(i,j,blockSizeX, blockSizeY, similarityGraph);
            }
        }
    }
    
    
    for(int i = 0; i<blockSizeX; ++i){
        for(int j = 0; j<blockSizeY; ++j){
            for(int k = 0; k<8; ++k){
                if(similarityGraph->at(i*blockSizeY+j)[k]){
                    int endX = i;
                    int endY = j;
                    getNeighbor(i, j, &endX, &endY, k);
                    if (similarityGraph->at(i*blockSizeY + j)[7] || similarityGraph->at(i*blockSizeY + j)[2] ||similarityGraph->at(i*blockSizeY + j)[0] || similarityGraph->at(i*blockSizeY + j)[5])
                        
                    
                    //reshapePixels(*similarityGraph, *rawGraphTest, blockSizeX, blockSizeY, pixelSize);
                     reshapePixelsGl(*similarityGraph, currentImageGraph, blockSizeX, blockSizeY, pixelSize);
                    //drawEdge(x, y, endX, endY, pixelSize, *rawGraphTest);
                }
            }
        }
        
//    }
//    
//    for(int i = 0; i<blockSizeX; ++i){
//        for(int j = 0; j<blockSizeY; ++j){
//            for(int k = 0; k<8; ++k){
//                if(similarityGraph->at(i*blockSizeY+j)[k]){
//                    int endX = i;
//                    int endY = j;
//                    getNeighbor(i, j, &endX, &endY, k);
//                    //reshapePixels(similarityGraph, *rawGraphTest, blockSizeX, blockSizeY, pixelSize);
//                    drawEdge(i, j, endX, endY, pixelSize, *rawGraphTest);
//                }
//            }
//        }
//        
    }
    
    return rawGraphTest;
}



/*
 * round each pixel to the nearest value in the new number of bits
 */
Image* ip_quantize_simple (Image* src, int bitsPerChannel)
{
    
    cerr << "This function is not implemented." << endl;
    return NULL;
}


/*
 * dither each pixel to the nearest value in the new number of bits
 * using a static 4x4 matrix
 */
Image* ip_quantize_ordered (Image* src, int bitsPerChannel)
{
    cerr << "This function is not implemented." << endl;
    return NULL;
}


/*
 * dither each pixel to the nearest value in the new number of bits
 * using error diffusion
 */
Image* ip_quantize_fs (Image* src, int bitsPerChannel)
{
    cerr << "This function is not implemented." << endl;
    return NULL;
}

/*
 * nearest neighbor sample
 */
Pixel ip_resample_nearest(Image* src, double x, double y) {
    cerr << "This function is not implemented." << endl;
    Pixel myPixel(0,0,0);
    
    return myPixel;
}

/*
 * bilinear sample
 */

Pixel ip_resample_bilinear(Image* src, double x, double y) {
    cerr << "This function is not implemented." << endl;
    Pixel myPixel(0,0,0);
    return myPixel;
}

/*
 * gausian sample
 */
Pixel ip_resample_gaussian(Image* src, double x, double y, int size, double sigma)
{
    cerr << "This function is not implemented." << endl;
    Pixel myPixel(0,0,0);
    return myPixel;
}

/*
 * rotate image using one of three sampling techniques
 */
Image* ip_rotate (Image* src, double theta, int x, int y, int samplingMode,
                  int gaussianFilterSize, double gaussianSigma)
{
    cerr << "This function is not implemented." << endl;
    return NULL;
}


/*
 * change saturation
 */
Image* ip_saturate (Image* src, double alpha)
{
    cerr << "This function is not implemented." << endl;
    return NULL;
}


/*
 * scale image using one of three sampling techniques
 */
Image* ip_scale (Image* src, double xFac, double yFac, int samplingMode,
                 int gaussianFilterSize, double gaussianSigma)
{
    cerr << "This function is not implemented." << endl;
    return NULL;
}


/*
 * threshold image
 */
Image* ip_threshold (Image* src, double cutoff)
{
    cerr << "This function is not implemented." << endl;
    return NULL;
}




