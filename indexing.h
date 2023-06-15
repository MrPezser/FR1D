//
// Created by Tsail on 6/14/2023.
//

#ifndef FR1D_INDEXING_H
#define FR1D_INDEXING_H

//indexing into state variable array u of degree nj for the jth node on the ith element
//column major i guess
#define iup(i, j, ndegr)  (((i)*(ndegr)) + (j))

#endif //FR1D_INDEXING_H
