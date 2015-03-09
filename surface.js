/* The files is assembled from the surface generating code of 3Dmol.js
 * This code was descended from GLmol, which implemented EDTSurf.
 * Acknowledgements follow. 
 * 
 * 
 * ProteinSurface.js by biochem_fan

Ported and modified for Javascript based on EDTSurf,
  whose license is as follows.

Permission to use, copy, modify, and distribute this program for any
purpose, with or without fee, is hereby granted, provided that this
copyright notice and the reference information appear in all copies or
substantial portions of the Software. It is provided "as is" without
express or implied warranty. 

Reference:
http://zhanglab.ccmb.med.umich.edu/EDTSurf/
D. Xu, Y. Zhang (2009) Generating Triangulated Macromolecular Surfaces
by Euclidean Distance Transform. PLoS ONE 4(12): e8140.

=======
 */

/*
 * These tables are based off those by Paul Bourke and Geoffrey Heller:
 * http://paulbourke.net/geometry/polygonise/
 * http://paulbourke.net/geometry/polygonise/table2.txt
 * 
 * However, they have been substantially modified to reflect a more 
 * sensible corner numbering scheme and the discrete nature of our voxel data
 * (resulting in fewer faces).
 */
var MarchingCube = MarchingCube || {};


//Encapsulate marching cube algorithm for isosurface generation
//(currently used by protein surface rendering and generic volumetric data reading)
MarchingCube = (function() {
 
 //Marching cube algorithm - assume data has been pre-treated so isovalue is 0 
 // (i.e. select points greater than 0)
 //origin -  vector of origin of volumetric data (default is (0,0,0))
 // nX, nY, nZ - specifies number of voxels in each dimension
 // scale - cube diagonal unit vector scale (3Dmol vector) (specifying distance between data points); diagonal of cube
 // - default is 1 - assumes unit cube (1,1,1) diag)
 // fulltable - if true, use full marching cubes and tritables - else use trimmed table (e.g. surf render)
 // voxel - if true, draws with a blocky voxel style (default false)
 // verts, faces - vertex and face arrays to fill up
 
 //to match with protein surface...
 var ISDONE = 2;
 var my = {};
 
 my.march = function(data, verts, faces, spec) {

     var fulltable = !!(spec.fulltable);
     var origin = (spec.hasOwnProperty('origin') && spec.origin.hasOwnProperty('x')) ? spec.origin : {x:0, y:0, z:0};
     var voxel = !!(spec.voxel);
     
     var nX = spec.nX || 0;
     var nY = spec.nY || 0;
     var nZ = spec.nZ || 0;
     
     var scale = spec.scale || 1.0;
     
     var unitCube = new $3Dmol.Vector3(1,1,1).multiplyScalar(scale);
     
     //keep track of calculated vertices to avoid repeats
     var vertnums = new Int32Array(nX*nY*nZ);
     
     var i, il;
     
     for (i = 0, il = vertnums.length; i < il; ++i)
         vertnums[i] = -1;

     // create (or retrieve) a vertex at the appropriate point for
     // the edge (p1,p2)
     
     var getVertex = function(i, j, k, code, p1, p2) {
         var pt = new $3Dmol.Vector3();
         pt.copy(origin);
         var val1 = !!(code & (1 << p1));
         var val2 = !!(code & (1 << p2));
          
         // p1 if they are the same or if !val1
         var p = p1;
         if (!val1 && val2)
             p = p2;
         
         // adjust i,j,k by p
         if (p & 1)
             k++;
         if (p & 2)
             j++;
         if (p & 4)
             i++;
 
         pt.x += unitCube.x*i;
         pt.y += unitCube.y*j;
         pt.z += unitCube.z*k;
 
         var index = ((nY * i) + j) * nZ + k;
         
         //Have to add option to do voxels
         if (!voxel) {
         
             if (vertnums[index] < 0) // not created yet
             {
                 vertnums[index] = verts.length;
                 verts.push( pt );
             }
             return vertnums[index];
         
         }
         
         else {
             verts.push(pt);
             return verts.length - 1;
         }
         
     };
         
     var intersects = new Int32Array(12);
     
     var etable = (fulltable) ? edgeTable2 : edgeTable;
     var tritable = (fulltable) ? triTable2 : triTable;
             
     //Run marching cubes algorithm
     for (i = 0; i < nX-1; ++i) {
         
         for (var j = 0; j < nY-1; ++j){
             
             for (var k = 0; k < nZ-1; ++k){
                 
                 var code = 0;
                 
                 for (var p = 0; p < 8; ++p) {
                     var index = ((nY * (i + ((p & 4) >> 2))) + j + ((p & 2) >> 1)) *
                                     nZ + k + (p & 1);

                     //TODO: Need to fix vpBits in protein surface for this to work
                     var val = !!(data[index] & ISDONE);
                     //var val = !!(data[index] > 0);   
                     
                     code |= val << p;                        
                 }
                 
                 if (code === 0 || code === 255)
                     continue;
                 
                 var ecode = etable[code];
                 
                 if (ecode === 0)
                     continue;
                     
                 var ttable = tritable[code];                        
                 
                 if (ecode & 1)
                     intersects[0] = getVertex(i, j, k, code, 0, 1);
                 if (ecode & 2)
                     intersects[1] = getVertex(i, j, k, code, 1, 3);
                 if (ecode & 4)
                     intersects[2] = getVertex(i, j, k, code, 3, 2);
                 if (ecode & 8)
                     intersects[3] = getVertex(i, j, k, code, 2, 0);
                 if (ecode & 16)
                     intersects[4] = getVertex(i, j, k, code, 4, 5);
                 if (ecode & 32)
                     intersects[5] = getVertex(i, j, k, code, 5, 7);
                 if (ecode & 64)
                     intersects[6] = getVertex(i, j, k, code, 7, 6);
                 if (ecode & 128)
                     intersects[7] = getVertex(i, j, k, code, 6, 4);
                 if (ecode & 256)
                     intersects[8] = getVertex(i, j, k, code, 0, 4);
                 if (ecode & 512)
                     intersects[9] = getVertex(i, j, k, code, 1, 5);
                 if (ecode & 1024)
                     intersects[10] = getVertex(i, j, k, code, 3, 7);
                 if (ecode & 2048)
                     intersects[11] = getVertex(i, j, k, code, 2, 6);       
                     
                 for (var t = 0; t < ttable.length; t += 3) {
                     
                     var a = intersects[ttable[t]],
                         b = intersects[ttable[t+1]],
                         c = intersects[ttable[t+2]];         
                                        
                     if (voxel && t >= 3) {
                         verts.push(verts[a]); a = verts.length - 1;
                         verts.push(verts[b]); b = verts.length - 1;
                         verts.push(verts[c]); c = verts.length - 1;
                     }

                     
                     faces.push(a); faces.push(b); faces.push(c);                               
                 }              
                 
             }
             
         }
         
     }
          
     
 };

 my.laplacianSmooth = function(numiter, verts, faces) {
         var tps = new Array(verts.length);
         var i, il, j, jl, k, kl;
         for (i = 0, il = verts.length; i < il; i++)
                 tps[i] = {
                     x : 0,
                     y : 0,
                     z : 0
                 };
         var vertdeg = new Array(20);
         var flagvert;
         for (i = 0; i < 20; i++)
                 vertdeg[i] = new Array(verts.length);
         for (i = 0, il = verts.length; i < il; i++)
                 vertdeg[0][i] = 0;
         for (i = 0, il = faces.length / 3; i < il; i++) {
             var aoffset = i*3, boffset = i*3 + 1, coffset = i*3 + 2;
             flagvert = true;
             for (j = 0, jl = vertdeg[0][faces[aoffset]]; j < jl; j++) {
                 if (faces[boffset] == vertdeg[j + 1][faces[aoffset]]) {
                     flagvert = false;
                     break;
                 }
             }
             if (flagvert) {
                 vertdeg[0][faces[aoffset]]++;
                 vertdeg[vertdeg[0][faces[aoffset]]][faces[aoffset]] = faces[boffset];
             }
             flagvert = true;
             for (j = 0, jl = vertdeg[0][faces[aoffset]]; j < jl; j++) {
                 if (faces[coffset] == vertdeg[j + 1][faces[aoffset]]) {
                     flagvert = false;
                     break;
                 }
             }
             if (flagvert) {
                 vertdeg[0][faces[aoffset]]++;
                 vertdeg[vertdeg[0][faces[aoffset]]][faces[aoffset]] = faces[coffset];
             }
             // b
             flagvert = true;
             for (j = 0, jl = vertdeg[0][faces[boffset]]; j < jl; j++) {
                 if (faces[aoffset] == vertdeg[j + 1][faces[boffset]]) {
                     flagvert = false;
                     break;
                 }
             }
             if (flagvert) {
                 vertdeg[0][faces[boffset]]++;
                 vertdeg[vertdeg[0][faces[boffset]]][faces[boffset]] = faces[aoffset];
             }
             flagvert = true;
             for (j = 0, jl = vertdeg[0][faces[boffset]]; j < jl; j++) {
                 if (faces[coffset] == vertdeg[j + 1][faces[boffset]]) {
                     flagvert = false;
                     break;
                 }
             }
             if (flagvert) {
                 vertdeg[0][faces[boffset]]++;
                 vertdeg[vertdeg[0][faces[boffset]]][faces[boffset]] = faces[coffset];
             }
             // c
             flagvert = true;
             for (j = 0; j < vertdeg[0][faces[coffset]]; j++) {
                 if (faces[aoffset] == vertdeg[j + 1][faces[coffset]]) {
                     flagvert = false;
                     break;
                 }
             }
             if (flagvert) {
                 vertdeg[0][faces[coffset]]++;
                 vertdeg[vertdeg[0][faces[coffset]]][faces[coffset]] = faces[aoffset];
             }
             flagvert = true;
             for (j = 0, jl = vertdeg[0][faces[coffset]]; j < jl; j++) {
                 if (faces[boffset] == vertdeg[j + 1][faces[coffset]]) {
                     flagvert = false;
                     break;
                 }
             }
             if (flagvert) {
                 vertdeg[0][faces[coffset]]++;
                 vertdeg[vertdeg[0][faces[coffset]]][faces[coffset]] = faces[boffset];
             }
         }

         var wt = 1.00;
         var wt2 = 0.50;
         var ssign;
         var scaleFactor = 1;
         var outwt = 0.75 / (scaleFactor + 3.5); // area-preserving
         for (k = 0; k < numiter; k++) {
                 for (i = 0, il = verts.length; i < il; i++) {
                         if (vertdeg[0][i] < 3) {
                                 tps[i].x = verts[i].x;
                                 tps[i].y = verts[i].y;
                                 tps[i].z = verts[i].z;
                         } else if (vertdeg[0][i] == 3 || vertdeg[0][i] == 4) {
                                 tps[i].x = 0;
                                 tps[i].y = 0;
                                 tps[i].z = 0;
                                 for (j = 0, jl = vertdeg[0][i]; j < jl; j++) {
                                         tps[i].x += verts[vertdeg[j + 1][i]].x;
                                         tps[i].y += verts[vertdeg[j + 1][i]].y;
                                         tps[i].z += verts[vertdeg[j + 1][i]].z;
                                 }
                                 tps[i].x += wt2 * verts[i].x;
                                 tps[i].y += wt2 * verts[i].y;
                                 tps[i].z += wt2 * verts[i].z;
                                 tps[i].x /= wt2 + vertdeg[0][i];
                                 tps[i].y /= wt2 + vertdeg[0][i];
                                 tps[i].z /= wt2 + vertdeg[0][i];
                         } else {
                                 tps[i].x = 0;
                                 tps[i].y = 0;
                                 tps[i].z = 0;
                                 for (j = 0, jl = vertdeg[0][i]; j < jl; j++) {
                                         tps[i].x += verts[vertdeg[j + 1][i]].x;
                                         tps[i].y += verts[vertdeg[j + 1][i]].y;
                                         tps[i].z += verts[vertdeg[j + 1][i]].z;
                                 }
                                 tps[i].x += wt * verts[i].x;
                                 tps[i].y += wt * verts[i].y;
                                 tps[i].z += wt * verts[i].z;
                                 tps[i].x /= wt + vertdeg[0][i];
                                 tps[i].y /= wt + vertdeg[0][i];
                                 tps[i].z /= wt + vertdeg[0][i];
                         }
                 }
                 for (i = 0, il = verts.length; i < il; i++) {
                         verts[i].x = tps[i].x;
                         verts[i].y = tps[i].y;
                         verts[i].z = tps[i].z;
                 }
                 /*
                  * computenorm(); for (var i = 0; i < vertnumber; i++) { if
                  * (verts[i].inout) ssign = 1; else ssign = -1; verts[i].x += ssign *
                  * outwt * verts[i].pn.x; verts[i].y += ssign * outwt *
                  * verts[i].pn.y; verts[i].z += ssign * outwt * verts[i].pn.z; }
                  */
         }
 };


 /*
  * These tables are based off those by Paul Bourke and Geoffrey Heller:
  * http://paulbourke.net/geometry/polygonise/
  * http://paulbourke.net/geometry/polygonise/table2.txt
  * 
  * However, they have been substantially modified to reflect a more 
  * sensible corner numbering scheme and the discrete nature of our voxel data
  * (resulting in fewer faces).
  */
 my.edgeTable = [ 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0,
         0xb00, 0x0, 0x0, 0x0, 0x700, 0x0, 0xd00, 0xe00, 0xf00, 0x0, 0x0, 0x0,
         0x8a, 0x0, 0x15, 0x0, 0x86, 0x0, 0x0, 0x0, 0x28c, 0x0, 0x813, 0xf19,
         0xe10, 0x0, 0x0, 0x0, 0x2a, 0x0, 0x0, 0x0, 0x126, 0x0, 0x0, 0x15, 0x1c,
         0x0, 0xf23, 0x419, 0xd20, 0x0, 0xa8, 0xa2, 0xaa, 0x0, 0x285, 0x9ab,
         0x8a2, 0x0, 0x2af, 0x125, 0xac, 0xfaa, 0xea3, 0xda9, 0xca0, 0x0, 0x0,
         0x0, 0x0, 0x0, 0x45, 0x0, 0x384, 0x0, 0x0, 0x0, 0x700, 0x8a, 0x83,
         0x648, 0x780, 0x0, 0x51, 0x0, 0x81a, 0x54, 0x55, 0x54, 0x56, 0x0, 0x51,
         0x0, 0xe5c, 0x14a, 0x451, 0x759, 0x650, 0x0, 0x0, 0x0, 0x2a, 0x0, 0x45,
         0x0, 0x1f6, 0x0, 0x0, 0x15, 0xdfc, 0x8a, 0x7f3, 0x4f9, 0x5f0, 0xb00,
         0x68, 0x921, 0x6a, 0x348, 0x245, 0x16f, 0x66, 0xb00, 0xe6f, 0xd65,
         0xc6c, 0x76a, 0x663, 0x569, 0x460, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0,
         0xf46, 0x0, 0x0, 0x45, 0x24c, 0x2a, 0x823, 0x29, 0xb40, 0x0, 0x0, 0x0,
         0x6ba, 0x0, 0x8f5, 0xfff, 0xef6, 0x0, 0xff, 0x2f5, 0x2fc, 0x9ea, 0x8f3,
         0xbf9, 0xaf0, 0x0, 0x0, 0x51, 0x152, 0x0, 0xf55, 0x45f, 0xd56, 0x54,
         0x357, 0x55, 0x154, 0x852, 0xb53, 0x59, 0x950, 0x700, 0x2c8, 0xc2,
         0x48a, 0xfc4, 0xec5, 0xdcf, 0xcc6, 0x2c4, 0x2cf, 0xc5, 0xcc, 0xbca,
         0xac3, 0x9c9, 0x8c0, 0x0, 0x0, 0x0, 0x0, 0xa8, 0x1a4, 0xa8, 0x7a6,
         0xa2, 0xa2, 0x2a4, 0xbac, 0xaa, 0xa3, 0x2a8, 0x3a0, 0xd00, 0xc18,
         0xd00, 0xe3a, 0x34, 0x35, 0x73f, 0x636, 0x924, 0x83f, 0xb35, 0xa3c,
         0x12a, 0x33, 0x339, 0x230, 0xe00, 0xe00, 0xc12, 0xd9a, 0x684, 0x795,
         0x49f, 0x596, 0x92, 0xb9f, 0x815, 0x99c, 0x9a, 0x393, 0x99, 0x190,
         0xf00, 0xe08, 0xd01, 0xc0a, 0x704, 0x605, 0x50f, 0x406, 0xb02, 0xa0f,
         0x905, 0x80c, 0x30a, 0x203, 0x109, 0x0 ];
 
 var edgeTable = new Uint32Array(my.edgeTable);
 
 var triTable = my.triTable = [ [], [], [], [], [], [], [], [ 11, 9, 8 ], [], [], [],
         [ 8, 10, 9 ], [], [ 10, 8, 11 ], [ 9, 11, 10 ],
         [ 8, 10, 9, 8, 11, 10 ], [], [], [], [ 1, 7, 3 ], [], [ 4, 2, 0 ], [],
         [ 2, 1, 7 ], [], [], [], [ 2, 7, 3, 2, 9, 7 ], [],
         [ 1, 4, 11, 1, 0, 4 ], [ 3, 8, 0, 11, 9, 4, 11, 10, 9 ],
         [ 4, 11, 9, 11, 10, 9 ], [], [], [], [ 5, 3, 1 ], [], [], [],
         [ 2, 5, 8, 2, 1, 5 ], [], [], [ 2, 4, 0 ], [ 3, 2, 4 ], [],
         [ 0, 9, 1, 8, 10, 5, 8, 11, 10 ], [ 3, 4, 0, 3, 10, 4 ],
         [ 5, 8, 10, 8, 11, 10 ], [], [ 3, 5, 7 ], [ 7, 1, 5 ],
         [ 1, 7, 3, 1, 5, 7 ], [], [ 9, 2, 0, 9, 7, 2 ],
         [ 0, 3, 8, 1, 7, 11, 1, 5, 7 ], [ 11, 1, 7, 1, 5, 7 ], [],
         [ 9, 1, 0, 5, 3, 2, 5, 7, 3 ], [ 8, 2, 5, 8, 0, 2 ],
         [ 2, 5, 3, 5, 7, 3 ], [ 3, 9, 1, 3, 8, 9, 7, 11, 10, 7, 10, 5 ],
         [ 9, 1, 0, 10, 7, 11, 10, 5, 7 ], [ 3, 8, 0, 7, 10, 5, 7, 11, 10 ],
         [ 11, 5, 7, 11, 10, 5 ], [], [], [], [], [], [ 0, 6, 2 ], [],
         [ 7, 2, 9, 7, 9, 8 ], [], [], [], [ 8, 10, 9 ], [ 7, 1, 3 ],
         [ 7, 1, 0 ], [ 6, 9, 3, 6, 10, 9 ], [ 7, 10, 8, 10, 9, 8 ], [],
         [ 6, 0, 4 ], [], [ 11, 1, 4, 11, 3, 1 ], [ 2, 4, 6 ],
         [ 2, 0, 4, 2, 4, 6 ], [ 2, 4, 6 ], [ 1, 4, 2, 4, 6, 2 ], [],
         [ 6, 0, 4 ], [], [ 2, 11, 3, 6, 9, 4, 6, 10, 9 ], [ 8, 6, 1, 8, 1, 3 ],
         [ 10, 0, 6, 0, 4, 6 ], [ 8, 0, 3, 9, 6, 10, 9, 4, 6 ],
         [ 10, 4, 6, 10, 9, 4 ], [], [], [], [ 5, 3, 1 ], [], [ 0, 6, 2 ], [],
         [ 7, 4, 8, 5, 2, 1, 5, 6, 2 ], [], [], [ 2, 4, 0 ],
         [ 7, 4, 8, 2, 11, 3, 10, 5, 6 ], [ 7, 1, 3 ],
         [ 5, 6, 10, 0, 9, 1, 8, 7, 4 ], [ 5, 6, 10, 7, 0, 3, 7, 4, 0 ],
         [ 10, 5, 6, 4, 8, 7 ], [ 9, 11, 8 ], [ 3, 5, 6 ],
         [ 0, 5, 11, 0, 11, 8 ], [ 6, 3, 5, 3, 1, 5 ], [ 3, 9, 6, 3, 8, 9 ],
         [ 9, 6, 0, 6, 2, 0 ], [ 0, 3, 8, 2, 5, 6, 2, 1, 5 ],
         [ 1, 6, 2, 1, 5, 6 ], [ 9, 11, 8 ], [ 1, 0, 9, 6, 10, 5, 11, 3, 2 ],
         [ 6, 10, 5, 2, 8, 0, 2, 11, 8 ], [ 3, 2, 11, 10, 5, 6 ],
         [ 10, 5, 6, 9, 3, 8, 9, 1, 3 ], [ 0, 9, 1, 5, 6, 10 ],
         [ 8, 0, 3, 10, 5, 6 ], [ 10, 5, 6 ], [], [], [], [], [], [], [],
         [ 1, 10, 2, 9, 11, 6, 9, 8, 11 ], [], [], [ 6, 0, 2 ],
         [ 3, 6, 9, 3, 2, 6 ], [ 3, 5, 1 ], [ 0, 5, 1, 0, 11, 5 ], [ 0, 3, 5 ],
         [ 6, 9, 11, 9, 8, 11 ], [], [], [], [ 4, 5, 9, 7, 1, 10, 7, 3, 1 ], [],
         [ 11, 6, 7, 2, 4, 5, 2, 0, 4 ],
         [ 11, 6, 7, 8, 0, 3, 1, 10, 2, 9, 4, 5 ],
         [ 6, 7, 11, 1, 10, 2, 9, 4, 5 ], [],
         [ 4, 1, 0, 4, 5, 1, 6, 7, 3, 6, 3, 2 ], [ 9, 4, 5, 0, 6, 7, 0, 2, 6 ],
         [ 4, 5, 9, 6, 3, 2, 6, 7, 3 ], [ 6, 7, 11, 5, 3, 8, 5, 1, 3 ],
         [ 6, 7, 11, 4, 1, 0, 4, 5, 1 ], [ 4, 5, 9, 3, 8, 0, 11, 6, 7 ],
         [ 9, 4, 5, 7, 11, 6 ], [], [], [ 0, 6, 4 ], [ 8, 6, 4, 8, 1, 6 ], [],
         [ 0, 10, 2, 0, 9, 10, 4, 8, 11, 4, 11, 6 ],
         [ 10, 2, 1, 6, 0, 3, 6, 4, 0 ], [ 10, 2, 1, 11, 4, 8, 11, 6, 4 ],
         [ 4, 2, 6 ], [ 1, 0, 9, 2, 4, 8, 2, 6, 4 ], [ 2, 4, 0, 2, 6, 4 ],
         [ 8, 2, 4, 2, 6, 4 ], [ 11, 4, 1, 11, 6, 4 ],
         [ 0, 9, 1, 4, 11, 6, 4, 8, 11 ], [ 3, 6, 0, 6, 4, 0 ],
         [ 8, 6, 4, 8, 11, 6 ], [ 10, 8, 9 ], [ 6, 3, 9, 6, 7, 3 ], [ 6, 7, 1 ],
         [ 10, 7, 1, 7, 3, 1 ], [ 7, 11, 6, 8, 10, 2, 8, 9, 10 ],
         [ 11, 6, 7, 10, 0, 9, 10, 2, 0 ], [ 2, 1, 10, 7, 11, 6, 8, 0, 3 ],
         [ 1, 10, 2, 6, 7, 11 ], [ 7, 2, 6, 7, 9, 2 ],
         [ 1, 0, 9, 3, 6, 7, 3, 2, 6 ], [ 7, 0, 6, 0, 2, 6 ],
         [ 2, 7, 3, 2, 6, 7 ], [ 7, 11, 6, 3, 9, 1, 3, 8, 9 ],
         [ 9, 1, 0, 11, 6, 7 ], [ 0, 3, 8, 11, 6, 7 ], [ 11, 6, 7 ], [], [], [],
         [], [ 5, 3, 7 ], [ 8, 5, 2, 8, 7, 5 ], [ 5, 3, 7 ],
         [ 1, 10, 2, 5, 8, 7, 5, 9, 8 ], [ 1, 7, 5 ], [ 1, 7, 5 ],
         [ 9, 2, 7, 9, 7, 5 ], [ 11, 3, 2, 8, 5, 9, 8, 7, 5 ],
         [ 1, 3, 7, 1, 7, 5 ], [ 0, 7, 1, 7, 5, 1 ], [ 9, 3, 5, 3, 7, 5 ],
         [ 9, 7, 5, 9, 8, 7 ], [ 8, 10, 11 ], [ 3, 4, 10, 3, 10, 11 ],
         [ 8, 10, 11 ], [ 5, 9, 4, 1, 11, 3, 1, 10, 11 ], [ 2, 4, 5 ],
         [ 5, 2, 4, 2, 0, 4 ], [ 0, 3, 8, 5, 9, 4, 10, 2, 1 ],
         [ 2, 1, 10, 9, 4, 5 ], [ 2, 8, 5, 2, 11, 8 ],
         [ 3, 2, 11, 1, 4, 5, 1, 0, 4 ], [ 9, 4, 5, 8, 2, 11, 8, 0, 2 ],
         [ 11, 3, 2, 9, 4, 5 ], [ 8, 5, 3, 5, 1, 3 ], [ 5, 0, 4, 5, 1, 0 ],
         [ 3, 8, 0, 4, 5, 9 ], [ 9, 4, 5 ], [ 11, 9, 10 ], [ 11, 9, 10 ],
         [ 1, 11, 4, 1, 10, 11 ], [ 8, 7, 4, 11, 1, 10, 11, 3, 1 ],
         [ 2, 7, 9, 2, 9, 10 ], [ 4, 8, 7, 0, 10, 2, 0, 9, 10 ],
         [ 2, 1, 10, 0, 7, 4, 0, 3, 7 ], [ 10, 2, 1, 8, 7, 4 ], [ 1, 7, 4 ],
         [ 3, 2, 11, 4, 8, 7, 9, 1, 0 ], [ 11, 4, 2, 4, 0, 2 ],
         [ 2, 11, 3, 7, 4, 8 ], [ 4, 1, 7, 1, 3, 7 ], [ 1, 0, 9, 8, 7, 4 ],
         [ 3, 4, 0, 3, 7, 4 ], [ 8, 7, 4 ], [ 8, 9, 10, 8, 10, 11 ],
         [ 3, 9, 11, 9, 10, 11 ], [ 0, 10, 8, 10, 11, 8 ],
         [ 10, 3, 1, 10, 11, 3 ], [ 2, 8, 10, 8, 9, 10 ], [ 9, 2, 0, 9, 10, 2 ],
         [ 8, 0, 3, 1, 10, 2 ], [ 10, 2, 1 ], [ 1, 11, 9, 11, 8, 9 ],
         [ 11, 3, 2, 0, 9, 1 ], [ 11, 0, 2, 11, 8, 0 ], [ 11, 3, 2 ],
         [ 8, 1, 3, 8, 9, 1 ], [ 9, 1, 0 ], [ 8, 0, 3 ], [] ];
  
 var edgeTable2 = [ 0x0, 0x109, 0x203, 0x30a, 0x80c, 0x905, 0xa0f,
         0xb06, 0x406, 0x50f, 0x605, 0x70c, 0xc0a, 0xd03, 0xe09, 0xf00, 0x190,
         0x99, 0x393, 0x29a, 0x99c, 0x895, 0xb9f, 0xa96, 0x596, 0x49f, 0x795,
         0x69c, 0xd9a, 0xc93, 0xf99, 0xe90, 0x230, 0x339, 0x33, 0x13a, 0xa3c,
         0xb35, 0x83f, 0x936, 0x636, 0x73f, 0x435, 0x53c, 0xe3a, 0xf33, 0xc39,
         0xd30, 0x3a0, 0x2a9, 0x1a3, 0xaa, 0xbac, 0xaa5, 0x9af, 0x8a6, 0x7a6,
         0x6af, 0x5a5, 0x4ac, 0xfaa, 0xea3, 0xda9, 0xca0, 0x8c0, 0x9c9, 0xac3,
         0xbca, 0xcc, 0x1c5, 0x2cf, 0x3c6, 0xcc6, 0xdcf, 0xec5, 0xfcc, 0x4ca,
         0x5c3, 0x6c9, 0x7c0, 0x950, 0x859, 0xb53, 0xa5a, 0x15c, 0x55, 0x35f,
         0x256, 0xd56, 0xc5f, 0xf55, 0xe5c, 0x55a, 0x453, 0x759, 0x650, 0xaf0,
         0xbf9, 0x8f3, 0x9fa, 0x2fc, 0x3f5, 0xff, 0x1f6, 0xef6, 0xfff, 0xcf5,
         0xdfc, 0x6fa, 0x7f3, 0x4f9, 0x5f0, 0xb60, 0xa69, 0x963, 0x86a, 0x36c,
         0x265, 0x16f, 0x66, 0xf66, 0xe6f, 0xd65, 0xc6c, 0x76a, 0x663, 0x569,
         0x460, 0x460, 0x569, 0x663, 0x76a, 0xc6c, 0xd65, 0xe6f, 0xf66, 0x66,
         0x16f, 0x265, 0x36c, 0x86a, 0x963, 0xa69, 0xb60, 0x5f0, 0x4f9, 0x7f3,
         0x6fa, 0xdfc, 0xcf5, 0xfff, 0xef6, 0x1f6, 0xff, 0x3f5, 0x2fc, 0x9fa,
         0x8f3, 0xbf9, 0xaf0, 0x650, 0x759, 0x453, 0x55a, 0xe5c, 0xf55, 0xc5f,
         0xd56, 0x256, 0x35f, 0x55, 0x15c, 0xa5a, 0xb53, 0x859, 0x950, 0x7c0,
         0x6c9, 0x5c3, 0x4ca, 0xfcc, 0xec5, 0xdcf, 0xcc6, 0x3c6, 0x2cf, 0x1c5,
         0xcc, 0xbca, 0xac3, 0x9c9, 0x8c0, 0xca0, 0xda9, 0xea3, 0xfaa, 0x4ac,
         0x5a5, 0x6af, 0x7a6, 0x8a6, 0x9af, 0xaa5, 0xbac, 0xaa, 0x1a3, 0x2a9,
         0x3a0, 0xd30, 0xc39, 0xf33, 0xe3a, 0x53c, 0x435, 0x73f, 0x636, 0x936,
         0x83f, 0xb35, 0xa3c, 0x13a, 0x33, 0x339, 0x230, 0xe90, 0xf99, 0xc93,
         0xd9a, 0x69c, 0x795, 0x49f, 0x596, 0xa96, 0xb9f, 0x895, 0x99c, 0x29a,
         0x393, 0x99, 0x190, 0xf00, 0xe09, 0xd03, 0xc0a, 0x70c, 0x605, 0x50f,
         0x406, 0xb06, 0xa0f, 0x905, 0x80c, 0x30a, 0x203, 0x109, 0x0 ];
  
 var triTable2 = [ [], [ 8, 3, 0 ], [ 9, 0, 1 ], [ 8, 3, 1, 8, 1, 9 ],
         [ 11, 2, 3 ], [ 11, 2, 0, 11, 0, 8 ], [ 11, 2, 3, 0, 1, 9 ],
         [ 2, 1, 11, 1, 9, 11, 11, 9, 8 ], [ 10, 1, 2 ], [ 8, 3, 0, 1, 2, 10 ],
         [ 9, 0, 2, 9, 2, 10 ], [ 3, 2, 8, 2, 10, 8, 8, 10, 9 ],
         [ 10, 1, 3, 10, 3, 11 ], [ 1, 0, 10, 0, 8, 10, 10, 8, 11 ],
         [ 0, 3, 9, 3, 11, 9, 9, 11, 10 ], [ 8, 10, 9, 8, 11, 10 ], [ 8, 4, 7 ],
         [ 3, 0, 4, 3, 4, 7 ], [ 1, 9, 0, 8, 4, 7 ],
         [ 9, 4, 1, 4, 7, 1, 1, 7, 3 ], [ 2, 3, 11, 7, 8, 4 ],
         [ 7, 11, 4, 11, 2, 4, 4, 2, 0 ], [ 3, 11, 2, 4, 7, 8, 9, 0, 1 ],
         [ 2, 7, 11, 2, 1, 7, 1, 4, 7, 1, 9, 4 ], [ 10, 1, 2, 8, 4, 7 ],
         [ 2, 10, 1, 0, 4, 7, 0, 7, 3 ], [ 4, 7, 8, 0, 2, 10, 0, 10, 9 ],
         [ 2, 7, 3, 2, 9, 7, 7, 9, 4, 2, 10, 9 ],
         [ 8, 4, 7, 11, 10, 1, 11, 1, 3 ],
         [ 11, 4, 7, 1, 4, 11, 1, 11, 10, 1, 0, 4 ],
         [ 3, 8, 0, 7, 11, 4, 11, 9, 4, 11, 10, 9 ],
         [ 7, 11, 4, 4, 11, 9, 11, 10, 9 ], [ 9, 5, 4 ], [ 3, 0, 8, 4, 9, 5 ],
         [ 5, 4, 0, 5, 0, 1 ], [ 4, 8, 5, 8, 3, 5, 5, 3, 1 ],
         [ 11, 2, 3, 9, 5, 4 ], [ 9, 5, 4, 8, 11, 2, 8, 2, 0 ],
         [ 3, 11, 2, 1, 5, 4, 1, 4, 0 ],
         [ 8, 5, 4, 2, 5, 8, 2, 8, 11, 2, 1, 5 ], [ 2, 10, 1, 9, 5, 4 ],
         [ 0, 8, 3, 5, 4, 9, 10, 1, 2 ], [ 10, 5, 2, 5, 4, 2, 2, 4, 0 ],
         [ 3, 4, 8, 3, 2, 4, 2, 5, 4, 2, 10, 5 ],
         [ 5, 4, 9, 1, 3, 11, 1, 11, 10 ],
         [ 0, 9, 1, 4, 8, 5, 8, 10, 5, 8, 11, 10 ],
         [ 3, 4, 0, 3, 10, 4, 4, 10, 5, 3, 11, 10 ],
         [ 4, 8, 5, 5, 8, 10, 8, 11, 10 ], [ 9, 5, 7, 9, 7, 8 ],
         [ 0, 9, 3, 9, 5, 3, 3, 5, 7 ], [ 8, 0, 7, 0, 1, 7, 7, 1, 5 ],
         [ 1, 7, 3, 1, 5, 7 ], [ 11, 2, 3, 8, 9, 5, 8, 5, 7 ],
         [ 9, 2, 0, 9, 7, 2, 2, 7, 11, 9, 5, 7 ],
         [ 0, 3, 8, 2, 1, 11, 1, 7, 11, 1, 5, 7 ],
         [ 2, 1, 11, 11, 1, 7, 1, 5, 7 ], [ 1, 2, 10, 5, 7, 8, 5, 8, 9 ],
         [ 9, 1, 0, 10, 5, 2, 5, 3, 2, 5, 7, 3 ],
         [ 5, 2, 10, 8, 2, 5, 8, 5, 7, 8, 0, 2 ],
         [ 10, 5, 2, 2, 5, 3, 5, 7, 3 ],
         [ 3, 9, 1, 3, 8, 9, 7, 11, 10, 7, 10, 5 ],
         [ 9, 1, 0, 10, 7, 11, 10, 5, 7 ], [ 3, 8, 0, 7, 10, 5, 7, 11, 10 ],
         [ 11, 5, 7, 11, 10, 5 ], [ 11, 7, 6 ], [ 0, 8, 3, 11, 7, 6 ],
         [ 9, 0, 1, 11, 7, 6 ], [ 7, 6, 11, 3, 1, 9, 3, 9, 8 ],
         [ 2, 3, 7, 2, 7, 6 ], [ 8, 7, 0, 7, 6, 0, 0, 6, 2 ],
         [ 1, 9, 0, 3, 7, 6, 3, 6, 2 ], [ 7, 6, 2, 7, 2, 9, 2, 1, 9, 7, 9, 8 ],
         [ 1, 2, 10, 6, 11, 7 ], [ 2, 10, 1, 7, 6, 11, 8, 3, 0 ],
         [ 11, 7, 6, 10, 9, 0, 10, 0, 2 ],
         [ 7, 6, 11, 3, 2, 8, 8, 2, 10, 8, 10, 9 ],
         [ 6, 10, 7, 10, 1, 7, 7, 1, 3 ],
         [ 6, 10, 1, 6, 1, 7, 7, 1, 0, 7, 0, 8 ],
         [ 9, 0, 3, 6, 9, 3, 6, 10, 9, 6, 3, 7 ],
         [ 6, 10, 7, 7, 10, 8, 10, 9, 8 ], [ 8, 4, 6, 8, 6, 11 ],
         [ 11, 3, 6, 3, 0, 6, 6, 0, 4 ], [ 0, 1, 9, 4, 6, 11, 4, 11, 8 ],
         [ 1, 9, 4, 11, 1, 4, 11, 3, 1, 11, 4, 6 ],
         [ 3, 8, 2, 8, 4, 2, 2, 4, 6 ], [ 2, 0, 4, 2, 4, 6 ],
         [ 1, 9, 0, 3, 8, 2, 2, 8, 4, 2, 4, 6 ], [ 9, 4, 1, 1, 4, 2, 4, 6, 2 ],
         [ 10, 1, 2, 11, 8, 4, 11, 4, 6 ],
         [ 10, 1, 2, 11, 3, 6, 6, 3, 0, 6, 0, 4 ],
         [ 0, 2, 10, 0, 10, 9, 4, 11, 8, 4, 6, 11 ],
         [ 2, 11, 3, 6, 9, 4, 6, 10, 9 ],
         [ 8, 4, 6, 8, 6, 1, 6, 10, 1, 8, 1, 3 ],
         [ 1, 0, 10, 10, 0, 6, 0, 4, 6 ], [ 8, 0, 3, 9, 6, 10, 9, 4, 6 ],
         [ 10, 4, 6, 10, 9, 4 ], [ 9, 5, 4, 7, 6, 11 ],
         [ 4, 9, 5, 3, 0, 8, 11, 7, 6 ], [ 6, 11, 7, 4, 0, 1, 4, 1, 5 ],
         [ 6, 11, 7, 4, 8, 5, 5, 8, 3, 5, 3, 1 ], [ 4, 9, 5, 6, 2, 3, 6, 3, 7 ],
         [ 9, 5, 4, 8, 7, 0, 0, 7, 6, 0, 6, 2 ],
         [ 4, 0, 1, 4, 1, 5, 6, 3, 7, 6, 2, 3 ], [ 7, 4, 8, 5, 2, 1, 5, 6, 2 ],
         [ 6, 11, 7, 1, 2, 10, 9, 5, 4 ],
         [ 11, 7, 6, 8, 3, 0, 1, 2, 10, 9, 5, 4 ],
         [ 11, 7, 6, 10, 5, 2, 2, 5, 4, 2, 4, 0 ],
         [ 7, 4, 8, 2, 11, 3, 10, 5, 6 ],
         [ 4, 9, 5, 6, 10, 7, 7, 10, 1, 7, 1, 3 ],
         [ 5, 6, 10, 0, 9, 1, 8, 7, 4 ], [ 5, 6, 10, 7, 0, 3, 7, 4, 0 ],
         [ 10, 5, 6, 4, 8, 7 ], [ 5, 6, 9, 6, 11, 9, 9, 11, 8 ],
         [ 0, 9, 5, 0, 5, 3, 3, 5, 6, 3, 6, 11 ],
         [ 0, 1, 5, 0, 5, 11, 5, 6, 11, 0, 11, 8 ],
         [ 11, 3, 6, 6, 3, 5, 3, 1, 5 ], [ 9, 5, 6, 3, 9, 6, 3, 8, 9, 3, 6, 2 ],
         [ 5, 6, 9, 9, 6, 0, 6, 2, 0 ], [ 0, 3, 8, 2, 5, 6, 2, 1, 5 ],
         [ 1, 6, 2, 1, 5, 6 ], [ 1, 2, 10, 5, 6, 9, 9, 6, 11, 9, 11, 8 ],
         [ 1, 0, 9, 6, 10, 5, 11, 3, 2 ], [ 6, 10, 5, 2, 8, 0, 2, 11, 8 ],
         [ 3, 2, 11, 10, 5, 6 ], [ 10, 5, 6, 9, 3, 8, 9, 1, 3 ],
         [ 0, 9, 1, 5, 6, 10 ], [ 8, 0, 3, 10, 5, 6 ], [ 10, 5, 6 ],
         [ 10, 6, 5 ], [ 8, 3, 0, 10, 6, 5 ], [ 0, 1, 9, 5, 10, 6 ],
         [ 10, 6, 5, 9, 8, 3, 9, 3, 1 ], [ 3, 11, 2, 10, 6, 5 ],
         [ 6, 5, 10, 2, 0, 8, 2, 8, 11 ], [ 1, 9, 0, 6, 5, 10, 11, 2, 3 ],
         [ 1, 10, 2, 5, 9, 6, 9, 11, 6, 9, 8, 11 ], [ 1, 2, 6, 1, 6, 5 ],
         [ 0, 8, 3, 2, 6, 5, 2, 5, 1 ], [ 5, 9, 6, 9, 0, 6, 6, 0, 2 ],
         [ 9, 6, 5, 3, 6, 9, 3, 9, 8, 3, 2, 6 ], [ 11, 6, 3, 6, 5, 3, 3, 5, 1 ],
         [ 0, 5, 1, 0, 11, 5, 5, 11, 6, 0, 8, 11 ],
         [ 0, 5, 9, 0, 3, 5, 3, 6, 5, 3, 11, 6 ],
         [ 5, 9, 6, 6, 9, 11, 9, 8, 11 ], [ 10, 6, 5, 4, 7, 8 ],
         [ 5, 10, 6, 7, 3, 0, 7, 0, 4 ], [ 5, 10, 6, 0, 1, 9, 8, 4, 7 ],
         [ 4, 5, 9, 6, 7, 10, 7, 1, 10, 7, 3, 1 ],
         [ 7, 8, 4, 2, 3, 11, 10, 6, 5 ],
         [ 11, 6, 7, 10, 2, 5, 2, 4, 5, 2, 0, 4 ],
         [ 11, 6, 7, 8, 0, 3, 1, 10, 2, 9, 4, 5 ],
         [ 6, 7, 11, 1, 10, 2, 9, 4, 5 ], [ 7, 8, 4, 5, 1, 2, 5, 2, 6 ],
         [ 4, 1, 0, 4, 5, 1, 6, 7, 3, 6, 3, 2 ],
         [ 9, 4, 5, 8, 0, 7, 0, 6, 7, 0, 2, 6 ], [ 4, 5, 9, 6, 3, 2, 6, 7, 3 ],
         [ 6, 7, 11, 4, 5, 8, 5, 3, 8, 5, 1, 3 ],
         [ 6, 7, 11, 4, 1, 0, 4, 5, 1 ], [ 4, 5, 9, 3, 8, 0, 11, 6, 7 ],
         [ 9, 4, 5, 7, 11, 6 ], [ 10, 6, 4, 10, 4, 9 ],
         [ 8, 3, 0, 9, 10, 6, 9, 6, 4 ], [ 1, 10, 0, 10, 6, 0, 0, 6, 4 ],
         [ 8, 6, 4, 8, 1, 6, 6, 1, 10, 8, 3, 1 ],
         [ 2, 3, 11, 6, 4, 9, 6, 9, 10 ],
         [ 0, 10, 2, 0, 9, 10, 4, 8, 11, 4, 11, 6 ],
         [ 10, 2, 1, 11, 6, 3, 6, 0, 3, 6, 4, 0 ],
         [ 10, 2, 1, 11, 4, 8, 11, 6, 4 ], [ 9, 1, 4, 1, 2, 4, 4, 2, 6 ],
         [ 1, 0, 9, 3, 2, 8, 2, 4, 8, 2, 6, 4 ], [ 2, 4, 0, 2, 6, 4 ],
         [ 3, 2, 8, 8, 2, 4, 2, 6, 4 ],
         [ 1, 4, 9, 11, 4, 1, 11, 1, 3, 11, 6, 4 ],
         [ 0, 9, 1, 4, 11, 6, 4, 8, 11 ], [ 11, 6, 3, 3, 6, 0, 6, 4, 0 ],
         [ 8, 6, 4, 8, 11, 6 ], [ 6, 7, 10, 7, 8, 10, 10, 8, 9 ],
         [ 9, 3, 0, 6, 3, 9, 6, 9, 10, 6, 7, 3 ],
         [ 6, 1, 10, 6, 7, 1, 7, 0, 1, 7, 8, 0 ],
         [ 6, 7, 10, 10, 7, 1, 7, 3, 1 ],
         [ 7, 11, 6, 3, 8, 2, 8, 10, 2, 8, 9, 10 ],
         [ 11, 6, 7, 10, 0, 9, 10, 2, 0 ], [ 2, 1, 10, 7, 11, 6, 8, 0, 3 ],
         [ 1, 10, 2, 6, 7, 11 ], [ 7, 2, 6, 7, 9, 2, 2, 9, 1, 7, 8, 9 ],
         [ 1, 0, 9, 3, 6, 7, 3, 2, 6 ], [ 8, 0, 7, 7, 0, 6, 0, 2, 6 ],
         [ 2, 7, 3, 2, 6, 7 ], [ 7, 11, 6, 3, 9, 1, 3, 8, 9 ],
         [ 9, 1, 0, 11, 6, 7 ], [ 0, 3, 8, 11, 6, 7 ], [ 11, 6, 7 ],
         [ 11, 7, 5, 11, 5, 10 ], [ 3, 0, 8, 7, 5, 10, 7, 10, 11 ],
         [ 9, 0, 1, 10, 11, 7, 10, 7, 5 ],
         [ 3, 1, 9, 3, 9, 8, 7, 10, 11, 7, 5, 10 ],
         [ 10, 2, 5, 2, 3, 5, 5, 3, 7 ],
         [ 5, 10, 2, 8, 5, 2, 8, 7, 5, 8, 2, 0 ],
         [ 9, 0, 1, 10, 2, 5, 5, 2, 3, 5, 3, 7 ],
         [ 1, 10, 2, 5, 8, 7, 5, 9, 8 ], [ 2, 11, 1, 11, 7, 1, 1, 7, 5 ],
         [ 0, 8, 3, 2, 11, 1, 1, 11, 7, 1, 7, 5 ],
         [ 9, 0, 2, 9, 2, 7, 2, 11, 7, 9, 7, 5 ],
         [ 11, 3, 2, 8, 5, 9, 8, 7, 5 ], [ 1, 3, 7, 1, 7, 5 ],
         [ 8, 7, 0, 0, 7, 1, 7, 5, 1 ], [ 0, 3, 9, 9, 3, 5, 3, 7, 5 ],
         [ 9, 7, 5, 9, 8, 7 ], [ 4, 5, 8, 5, 10, 8, 8, 10, 11 ],
         [ 3, 0, 4, 3, 4, 10, 4, 5, 10, 3, 10, 11 ],
         [ 0, 1, 9, 4, 5, 8, 8, 5, 10, 8, 10, 11 ],
         [ 5, 9, 4, 1, 11, 3, 1, 10, 11 ],
         [ 3, 8, 4, 3, 4, 2, 2, 4, 5, 2, 5, 10 ],
         [ 10, 2, 5, 5, 2, 4, 2, 0, 4 ], [ 0, 3, 8, 5, 9, 4, 10, 2, 1 ],
         [ 2, 1, 10, 9, 4, 5 ], [ 8, 4, 5, 2, 8, 5, 2, 11, 8, 2, 5, 1 ],
         [ 3, 2, 11, 1, 4, 5, 1, 0, 4 ], [ 9, 4, 5, 8, 2, 11, 8, 0, 2 ],
         [ 11, 3, 2, 9, 4, 5 ], [ 4, 5, 8, 8, 5, 3, 5, 1, 3 ],
         [ 5, 0, 4, 5, 1, 0 ], [ 3, 8, 0, 4, 5, 9 ], [ 9, 4, 5 ],
         [ 7, 4, 11, 4, 9, 11, 11, 9, 10 ],
         [ 3, 0, 8, 7, 4, 11, 11, 4, 9, 11, 9, 10 ],
         [ 11, 7, 4, 1, 11, 4, 1, 10, 11, 1, 4, 0 ],
         [ 8, 7, 4, 11, 1, 10, 11, 3, 1 ],
         [ 2, 3, 7, 2, 7, 9, 7, 4, 9, 2, 9, 10 ],
         [ 4, 8, 7, 0, 10, 2, 0, 9, 10 ], [ 2, 1, 10, 0, 7, 4, 0, 3, 7 ],
         [ 10, 2, 1, 8, 7, 4 ], [ 2, 11, 7, 2, 7, 1, 1, 7, 4, 1, 4, 9 ],
         [ 3, 2, 11, 4, 8, 7, 9, 1, 0 ], [ 7, 4, 11, 11, 4, 2, 4, 0, 2 ],
         [ 2, 11, 3, 7, 4, 8 ], [ 9, 1, 4, 4, 1, 7, 1, 3, 7 ],
         [ 1, 0, 9, 8, 7, 4 ], [ 3, 4, 0, 3, 7, 4 ], [ 8, 7, 4 ],
         [ 8, 9, 10, 8, 10, 11 ], [ 0, 9, 3, 3, 9, 11, 9, 10, 11 ],
         [ 1, 10, 0, 0, 10, 8, 10, 11, 8 ], [ 10, 3, 1, 10, 11, 3 ],
         [ 3, 8, 2, 2, 8, 10, 8, 9, 10 ], [ 9, 2, 0, 9, 10, 2 ],
         [ 8, 0, 3, 1, 10, 2 ], [ 10, 2, 1 ], [ 2, 11, 1, 1, 11, 9, 11, 8, 9 ],
         [ 11, 3, 2, 0, 9, 1 ], [ 11, 0, 2, 11, 8, 0 ], [ 11, 3, 2 ],
         [ 8, 1, 3, 8, 9, 1 ], [ 9, 1, 0 ], [ 8, 0, 3 ], [] ];
         
         return my;
})();




var ProteinSurface = function() {

    // constants for vpbits bitmasks
    /** @const */
    var INOUT = 1;
    /** @const */
    var ISDONE = 2;
    /** @const */
    var ISBOUND = 4;

    var ptranx = 0, ptrany = 0, ptranz = 0;
    var probeRadius = 1.4;
    var defaultScaleFactor = 2;
    var scaleFactor = defaultScaleFactor; // 2 is .5A grid; if this is made user configurable,
                            // also have to adjust offset used to find non-shown
                            // atoms
    var pHeight = 0, pWidth = 0, pLength = 0;
    var cutRadius = 0;
    var vpBits = null; // uint8 array of bitmasks
    var vpDistance = null; // floatarray of _squared_ distances
    var vpAtomID = null; // intarray
    var vertnumber = 0, facenumber = 0;
    var pminx = 0, pminy = 0, pminz = 0, pmaxx = 0, pmaxy = 0, pmaxz = 0;

    var vdwRadii = {
            "H" : 1.2,
            "Li" : 1.82,
            "Na" : 2.27,
            "K" : 2.75,
            "C" : 1.7,
            "N" : 1.55,
            "O" : 1.52,
            "F" : 1.47,
            "P" : 1.80,
            "S" : 1.80,
            "CL" : 1.75,
            "BR" : 1.85,
            "SE" : 1.90,
            "ZN" : 1.39,
            "CU" : 1.4,
            "NI" : 1.63,
            "X" : 2
        };
    
    /** @param {AtomSpec} atom */
    var getVDWIndex = function(atom) {
        if(!atom.elem || typeof(vdwRadii[atom.elem]) == "undefined") {
            return "X";
        }
        return atom.elem;
    };
    
    var depty = {}, widxz = {};
    var faces, verts;
    var nb = [ new Int32Array([ 1, 0, 0 ]), new Int32Array([ -1, 0, 0 ]), 
               new Int32Array([ 0, 1, 0 ]), new Int32Array([ 0, -1, 0 ]),
               new Int32Array([ 0, 0, 1 ]), 
               new Int32Array([ 0, 0, -1 ]), 
               new Int32Array([ 1, 1, 0 ]), 
               new Int32Array([ 1, -1, 0 ]), 
               new Int32Array([ -1, 1, 0 ]),
               new Int32Array([ -1, -1, 0 ]), 
               new Int32Array([ 1, 0, 1 ]), 
               new Int32Array([ 1, 0, -1 ]), 
               new Int32Array([ -1, 0, 1 ]),
               new Int32Array([ -1, 0, -1 ]), 
               new Int32Array([ 0, 1, 1 ]), 
               new Int32Array([ 0, 1, -1 ]), 
               new Int32Array([ 0, -1, 1 ]),
               new Int32Array([ 0, -1, -1 ]), 
               new Int32Array([ 1, 1, 1 ]), 
               new Int32Array([ 1, 1, -1 ]), 
               new Int32Array([ 1, -1, 1 ]),
               new Int32Array([ -1, 1, 1 ]), 
               new Int32Array([ 1, -1, -1 ]), 
               new Int32Array([ -1, -1, 1 ]), 
               new Int32Array([ -1, 1, -1 ]),
               new Int32Array([ -1, -1, -1 ]) ];

    var origextent;

    var inOrigExtent = function(x, y, z) {
        if (x < origextent[0][0] || x > origextent[1][0])
            return false;
        if (y < origextent[0][1] || y > origextent[1][1])
            return false;
        if (z < origextent[0][2] || z > origextent[1][2])
            return false;
        return true;
    };

    this.getFacesAndVertices = function(atomlist) {
        var atomsToShow = {};
        var i, il;
        for (i = 0, il = atomlist.length; i < il; i++)
            atomsToShow[atomlist[i]] = true;
        var vertices = verts;
        for (i = 0, il = vertices.length; i < il; i++) {
            vertices[i].x = vertices[i].x / scaleFactor - ptranx;
            vertices[i].y = vertices[i].y / scaleFactor - ptrany;
            vertices[i].z = vertices[i].z / scaleFactor - ptranz;
        }

        var finalfaces = [];
        for (i = 0, il = faces.length; i < il; i += 3) {
            //var f = faces[i];
            var fa = faces[i], fb = faces[i+1], fc = faces[i+2];
            var a = vertices[fa]['atomid'], b = vertices[fb]['atomid'], c = vertices[fc]['atomid'];

            // must be a unique face for each atom
            var which = a;
            if (b < which)
                which = b;
            if (c < which)
                which = c;
            if (!atomsToShow[which]) {
                continue;
            }
            var av = vertices[faces[i]];
            var bv = vertices[faces[i+1]];
            var cv = vertices[faces[i+2]];

            if (fa !== fb && fb !== fc && fa !== fc){
                finalfaces.push(fa); 
                finalfaces.push(fb); 
                finalfaces.push(fc); 
            }
               
        }

        //try to help the garbage collector
        vpBits = null; // uint8 array of bitmasks
        vpDistance = null; // floatarray
        vpAtomID = null; // intarray
        
        return {
            'vertices' : vertices,
            'faces' : finalfaces
        };
    };


    this.initparm = function(extent, btype, volume) {
        
        var margin = (1 / scaleFactor) * 5.5; // need margin to avoid
                                                // boundary/round off effects
        origextent = extent;
        pminx = extent[0][0]; pmaxx = extent[1][0];
        pminy = extent[0][1]; pmaxy = extent[1][1];
        pminz = extent[0][2]; pmaxz = extent[1][2];

        if (!btype) {
            pminx -= margin;
            pminy -= margin;
            pminz -= margin;
            pmaxx += margin;
            pmaxy += margin;
            pmaxz += margin;
        } else {
            pminx -= probeRadius + margin;
            pminy -= probeRadius + margin;
            pminz -= probeRadius + margin;
            pmaxx += probeRadius + margin;
            pmaxy += probeRadius + margin;
            pmaxz += probeRadius + margin;
        }

        pminx = Math.floor(pminx * scaleFactor) / scaleFactor;
        pminy = Math.floor(pminy * scaleFactor) / scaleFactor;
        pminz = Math.floor(pminz * scaleFactor) / scaleFactor;
        pmaxx = Math.ceil(pmaxx * scaleFactor) / scaleFactor;
        pmaxy = Math.ceil(pmaxy * scaleFactor) / scaleFactor;
        pmaxz = Math.ceil(pmaxz * scaleFactor) / scaleFactor;

        ptranx = -pminx;
        ptrany = -pminy;
        ptranz = -pminz;

        pLength = Math.ceil(scaleFactor * (pmaxx - pminx)) + 1;
        pWidth = Math.ceil(scaleFactor * (pmaxy - pminy)) + 1;
        pHeight = Math.ceil(scaleFactor * (pmaxz - pminz)) + 1;

        this.boundingatom(btype);
        cutRadius = probeRadius * scaleFactor;

        vpBits = new Uint8Array(pLength * pWidth * pHeight);
        vpDistance = new Float64Array(pLength * pWidth * pHeight); // float 32
        // doesn't
        // play
        // nicely
        // with
        // native
        // floats
        vpAtomID = new Int32Array(pLength * pWidth * pHeight);
        console.log("Box size: ", pLength, pWidth, pHeight, vpBits.length);
    };

    this.boundingatom = function(btype) {
        var tradius = [];
        var txz, tdept, sradius, idx;
        flagradius = btype;

        for ( var i in vdwRadii) {
            if(!vdwRadii.hasOwnProperty(i))
                continue;
            var r = vdwRadii[i];
            if (!btype)
                tradius[i] = r * scaleFactor + 0.5;
            else
                tradius[i] = (r + probeRadius) * scaleFactor + 0.5;

            sradius = tradius[i] * tradius[i];
            widxz[i] = Math.floor(tradius[i]) + 1;
            depty[i] = new Int32Array(widxz[i] * widxz[i]);
            indx = 0;
            for (j = 0; j < widxz[i]; j++) {
                for (k = 0; k < widxz[i]; k++) {
                    txz = j * j + k * k;
                    if (txz > sradius)
                        depty[i][indx] = -1; // outside
                    else {
                        tdept = Math.sqrt(sradius - txz);
                        depty[i][indx] = Math.floor(tdept);
                    }
                    indx++;
                }
            }
        }
    };

    this.fillvoxels = function(atoms, atomlist) { // (int seqinit,int
        // seqterm,bool
        // atomtype,atom*
        // proseq,bool bcolor)
        var i, il;
        for (i = 0, il = vpBits.length; i < il; i++) {
            vpBits[i] = 0;
            vpDistance[i] = -1.0;
            vpAtomID[i] = -1;
        }

        for (i in atomlist) {
            var atom = atoms[atomlist[i]];
            if (atom === undefined)
                continue;
            this.fillAtom(atom, atoms);
        }

        for (i = 0, il = vpBits.length; i < il; i++)
            if (vpBits[i] & INOUT)
                vpBits[i] |= ISDONE;

    };


    this.fillAtom = function(atom, atoms) {
        var cx, cy, cz, ox, oy, oz, mi, mj, mk, i, j, k, si, sj, sk;
        var ii, jj, kk, n;
        cx = Math.floor(0.5 + scaleFactor * (atom.x + ptranx));
        cy = Math.floor(0.5 + scaleFactor * (atom.y + ptrany));
        cz = Math.floor(0.5 + scaleFactor * (atom.z + ptranz));

        var at = getVDWIndex(atom);
        var nind = 0;
        var cnt = 0;
        var pWH = pWidth*pHeight;
        
        for (i = 0, n = widxz[at]; i < n; i++) {
            for (j = 0; j < n; j++) {
                if (depty[at][nind] != -1) {
                    for (ii = -1; ii < 2; ii++) {
                        for (jj = -1; jj < 2; jj++) {
                            for (kk = -1; kk < 2; kk++) {
                                if (ii !== 0 && jj !== 0 && kk !== 0) {
                                    mi = ii * i;
                                    mk = kk * j;
                                    for (k = 0; k <= depty[at][nind]; k++) {
                                        mj = k * jj;
                                        si = cx + mi;
                                        sj = cy + mj;
                                        sk = cz + mk;
                                        if (si < 0 || sj < 0 || 
                                                sk < 0 ||
                                                si >= pLength || 
                                                sj >= pWidth || 
                                                sk >= pHeight)
                                            continue;
                                        var index = si * pWH + sj * pHeight + sk;

                                        if (!(vpBits[index] & INOUT)) {
                                            vpBits[index] |= INOUT;
                                            vpAtomID[index] = atom.serial;
                                        } else {
                                            var atom2 = atoms[vpAtomID[index]];
                                            ox = Math.floor(0.5 + scaleFactor *
                                                    (atom2.x + ptranx));
                                            oy = Math.floor(0.5 + scaleFactor *
                                                    (atom2.y + ptrany));
                                            oz = Math.floor(0.5 + scaleFactor *
                                                    (atom2.z + ptranz));
                                            if (mi * mi + mj * mj + mk * mk < ox *
                                                    ox + oy * oy + oz * oz)
                                                vpAtomID[index] = atom.serial;
                                        }

                                    }// k
                                }// if
                            }// kk
                        }// jj
                    }// ii
                }// if
                nind++;
            }// j
        }// i
    };

    this.fillvoxelswaals = function(atoms, atomlist) {
        var i, il;
        for (i = 0, il = vpBits.length; i < il; i++)
            vpBits[i] &= ~ISDONE; // not isdone

        for (i in atomlist) {
            var atom = atoms[atomlist[i]];
            if (atom === undefined)
                continue;

            this.fillAtomWaals(atom, atoms);
        }
    };

    this.fillAtomWaals = function(atom, atoms) {
        var cx, cy, cz, ox, oy, oz, nind = 0;
        var mi, mj, mk, si, sj, sk, i, j, k, ii, jj, kk, n;
        cx = Math.floor(0.5 + scaleFactor * (atom.x + ptranx));
        cy = Math.floor(0.5 + scaleFactor * (atom.y + ptrany));
        cz = Math.floor(0.5 + scaleFactor * (atom.z + ptranz));

        var at = getVDWIndex(atom);
        var pWH = pWidth*pHeight;
        for (i = 0, n = widxz[at]; i < n; i++) {
            for (j = 0; j < n; j++) {
                if (depty[at][nind] != -1) {
                    for (ii = -1; ii < 2; ii++) {
                        for (jj = -1; jj < 2; jj++) {
                            for (kk = -1; kk < 2; kk++) {
                                if (ii !== 0 && jj !== 0 && kk !== 0) {
                                    mi = ii * i;
                                    mk = kk * j;
                                    for (k = 0; k <= depty[at][nind]; k++) {
                                        mj = k * jj;
                                        si = cx + mi;
                                        sj = cy + mj;
                                        sk = cz + mk;
                                        if (si < 0 || sj < 0 || 
                                                sk < 0 || 
                                                si >= pLength || 
                                                sj >= pWidth || 
                                                sk >= pHeight)
                                            continue;
                                        var index = si * pWH + sj * pHeight + sk;
                                        if (!(vpBits[index] & ISDONE)) {
                                            vpBits[index] |= ISDONE;
                                            vpAtomID[index] = atom.serial;
                                        } else {
                                            var atom2 = atoms[vpAtomID[index]];
                                            ox = Math.floor(0.5 + scaleFactor * (atom2.x + ptranx));
                                            oy = Math.floor(0.5 + scaleFactor * (atom2.y + ptrany));
                                            oz = Math.floor(0.5 + scaleFactor * (atom2.z + ptranz));
                                            if (mi * mi + mj * mj + mk * mk < ox *
                                                    ox + oy * oy + oz * oz)
                                                vpAtomID[index] = atom.serial;
                                        }

                                    }// k
                                }// if
                            }// kk
                        }// jj
                    }// ii
                }// if
                nind++;
            }// j
        }// i
    };

    this.buildboundary = function() {
        var pWH = pWidth*pHeight;
        for (i = 0; i < pLength; i++) {
            for (j = 0; j < pHeight; j++) {
                for (k = 0; k < pWidth; k++) {
                    var index = i * pWH + k * pHeight + j;
                    if (vpBits[index] & INOUT) {
                        var flagbound = false;
                        var ii = 0;
                        while (ii < 26) {
                            var ti = i + nb[ii][0], tj = j + nb[ii][2], tk = k +
                                    nb[ii][1];
                            if (ti > -1 && 
                                ti < pLength && 
                                tk > -1 && 
                                tk < pWidth && 
                                tj > -1 && 
                                tj < pHeight && 
                                !(vpBits[ti * pWH + tk * pHeight + tj] & INOUT)) {
                                vpBits[index] |= ISBOUND;
                                break;
                            } else
                                ii++;
                        }
                    }
                }
            }
        }
    };

    // a little class for 3d array, should really generalize this and
    // use throughout...
    var PointGrid = function(length, width, height) {
        // the standard says this is zero initialized
        var data = new Int32Array(length * width * height * 3);

        // set position x,y,z to pt, which has ix,iy,and iz
        this.set = function(x, y, z, pt) {
            var index = ((((x * width) + y) * height) + z) * 3;
            data[index] = pt.ix;
            data[index + 1] = pt.iy;
            data[index + 2] = pt.iz;
        };

        // return point at x,y,z
        this.get = function(x, y, z) {
            var index = ((((x * width) + y) * height) + z) * 3;
            return {
                ix : data[index],
                iy : data[index + 1],
                iz : data[index + 2]
            };
        };
    };

    this.fastdistancemap = function() {
        var eliminate = 0;
        var certificate;
        var i, j, k, n;

        var boundPoint = new PointGrid(pLength, pWidth, pHeight);
        var pWH = pWidth*pHeight;
        var cutRSq = cutRadius*cutRadius;
        
        var inarray = [];
        var outarray = [];
        
        var index;
        
        for (i = 0; i < pLength; i++) {
            for (j = 0; j < pWidth; j++) {
                for (k = 0; k < pHeight; k++) {
                    index = i * pWH + j * pHeight + k;
                    vpBits[index] &= ~ISDONE; // isdone = false
                    if (vpBits[index] & INOUT) {
                        if (vpBits[index] & ISBOUND) {
                            var triple = {
                                ix : i,
                                iy : j,
                                iz : k
                            };
                            boundPoint.set(i, j, k, triple);
                            inarray.push(triple);
                            vpDistance[index] = 0;
                            vpBits[index] |= ISDONE;
                            vpBits[index] &= ~ISBOUND;
                        } 
                    }
                }
            }
        }

        do {
            outarray = this.fastoneshell(inarray, boundPoint);
            inarray = [];
            for (i = 0, n = outarray.length; i < n; i++) {
                index = pWH * outarray[i].ix + pHeight *
                    outarray[i].iy + outarray[i].iz;
                vpBits[index] &= ~ISBOUND;
                if (vpDistance[index] <= 1.0404 * cutRSq) {
                    inarray.push({
                        ix : outarray[i].ix,
                        iy : outarray[i].iy,
                        iz : outarray[i].iz
                    });
                }
            }
        } while (inarray.length !== 0);

        inarray = [];
        outarray = [];
        boundPoint = null;
        
        var cutsf = scaleFactor - 0.5;
        if (cutsf < 0)
            cutsf = 0;
        var cutoff = cutRSq - 0.50 / (0.1 + cutsf);
        for (i = 0; i < pLength; i++) {
            for (j = 0; j < pWidth; j++) {
                for (k = 0; k < pHeight; k++) {
                    index = i * pWH + j * pHeight + k;
                    vpBits[index] &= ~ISBOUND;
                    // ses solid
                    if (vpBits[index] & INOUT) {
                        if (!(vpBits[index] & ISDONE) ||
                                ((vpBits[index] & ISDONE) && vpDistance[index] >= cutoff)) {
                            vpBits[index] |= ISBOUND;
                        }
                    }
                }
            }
        }

    };

    this.fastoneshell = function(inarray, boundPoint) { // (int* innum,int
        // *allocout,voxel2
        // ***boundPoint, int*
        // outnum, int *elimi)
        var tx, ty, tz;
        var dx, dy, dz;
        var i, j, n;
        var square;
        var bp, index;
        var outarray = [];
        if (inarray.length === 0)
            return outarray;

        tnv = {
            ix : -1,
            iy : -1,
            iz : -1
        };
        var pWH = pWidth*pHeight;
        for ( i = 0, n = inarray.length; i < n; i++) {
            tx = inarray[i].ix;
            ty = inarray[i].iy;
            tz = inarray[i].iz;
            bp = boundPoint.get(tx, ty, tz);

            for (j = 0; j < 6; j++) {
                tnv.ix = tx + nb[j][0];
                tnv.iy = ty + nb[j][1];
                tnv.iz = tz + nb[j][2];
                
                if (tnv.ix < pLength && tnv.ix > -1 && tnv.iy < pWidth &&
                        tnv.iy > -1 && tnv.iz < pHeight && tnv.iz > -1) {
                    index = tnv.ix * pWH + pHeight * tnv.iy + tnv.iz;
                    
                    if ((vpBits[index] & INOUT) && !(vpBits[index] & ISDONE)) {
    
                        boundPoint.set(tnv.ix, tnv.iy, tz + nb[j][2], bp);
                        dx = tnv.ix - bp.ix;
                        dy = tnv.iy - bp.iy;
                        dz = tnv.iz - bp.iz;
                        square = dx * dx + dy * dy + dz * dz;
                        vpDistance[index] = square;
                        vpBits[index] |= ISDONE;
                        vpBits[index] |= ISBOUND;
    
                        outarray.push({
                            ix : tnv.ix,
                            iy : tnv.iy,
                            iz : tnv.iz
                        });
                    } else if ((vpBits[index] & INOUT) && (vpBits[index] & ISDONE)) {
    
                        dx = tnv.ix - bp.ix;
                        dy = tnv.iy - bp.iy;
                        dz = tnv.iz - bp.iz;
                        square = dx * dx + dy * dy + dz * dz;
                        if (square < vpDistance[index]) {
                            boundPoint.set(tnv.ix, tnv.iy, tnv.iz, bp);
    
                            vpDistance[index] = square;
                            if (!(vpBits[index] & ISBOUND)) {
                                vpBits[index] |= ISBOUND;
                                outarray.push({
                                    ix : tnv.ix,
                                    iy : tnv.iy,
                                    iz : tnv.iz
                                });
                            }
                        }
                    }
                }
            }
        }

        // console.log("part1", positout);

        for (i = 0, n = inarray.length; i < n; i++) {
            tx = inarray[i].ix;
            ty = inarray[i].iy;
            tz = inarray[i].iz;
            bp = boundPoint.get(tx, ty, tz);

            for (j = 6; j < 18; j++) {
                tnv.ix = tx + nb[j][0];
                tnv.iy = ty + nb[j][1];
                tnv.iz = tz + nb[j][2];

                if(tnv.ix < pLength && tnv.ix > -1 && tnv.iy < pWidth &&
                        tnv.iy > -1 && tnv.iz < pHeight && tnv.iz > -1) {
                    index = tnv.ix * pWH + pHeight * tnv.iy + tnv.iz;
                    
                    if ((vpBits[index] & INOUT) && !(vpBits[index] & ISDONE)) {
                        boundPoint.set(tnv.ix, tnv.iy, tz + nb[j][2], bp);
    
                        dx = tnv.ix - bp.ix;
                        dy = tnv.iy - bp.iy;
                        dz = tnv.iz - bp.iz;
                        square = dx * dx + dy * dy + dz * dz;
                        vpDistance[index] = square;
                        vpBits[index] |= ISDONE;
                        vpBits[index] |= ISBOUND;
    
                        outarray.push({
                            ix : tnv.ix,
                            iy : tnv.iy,
                            iz : tnv.iz
                        });
                    } else if ((vpBits[index] & INOUT) && (vpBits[index] & ISDONE)) {
                        dx = tnv.ix - bp.ix;
                        dy = tnv.iy - bp.iy;
                        dz = tnv.iz - bp.iz;
                        square = dx * dx + dy * dy + dz * dz;
                        if (square < vpDistance[index]) {
                            boundPoint.set(tnv.ix, tnv.iy, tnv.iz, bp);
                            vpDistance[index] = square;
                            if (!(vpBits[index] & ISBOUND)) {
                                vpBits[index] |= ISBOUND;
                                outarray.push({
                                    ix : tnv.ix,
                                    iy : tnv.iy,
                                    iz : tnv.iz
                                });
                            }
                        }
                    }
                }
            }
        }

        // console.log("part2", positout);

        for (i = 0, n = inarray.length; i < n; i++) {
            tx = inarray[i].ix;
            ty = inarray[i].iy;
            tz = inarray[i].iz;
            bp = boundPoint.get(tx, ty, tz);

            for (j = 18; j < 26; j++) {
                tnv.ix = tx + nb[j][0];
                tnv.iy = ty + nb[j][1];
                tnv.iz = tz + nb[j][2];

                if (tnv.ix < pLength && tnv.ix > -1 && tnv.iy < pWidth &&
                        tnv.iy > -1 && tnv.iz < pHeight && tnv.iz > -1) {
                    index = tnv.ix * pWH + pHeight * tnv.iy + tnv.iz;

                    if ((vpBits[index] & INOUT) && !(vpBits[index] & ISDONE)) {
                        boundPoint.set(tnv.ix, tnv.iy, tz + nb[j][2], bp);

                        dx = tnv.ix - bp.ix;
                        dy = tnv.iy - bp.iy;
                        dz = tnv.iz - bp.iz;
                        square = dx * dx + dy * dy + dz * dz;
                        vpDistance[index] = square;
                        vpBits[index] |= ISDONE;
                        vpBits[index] |= ISBOUND;

                        outarray.push({
                            ix : tnv.ix,
                            iy : tnv.iy,
                            iz : tnv.iz
                        });
                    } else if ((vpBits[index] & INOUT)  && (vpBits[index] & ISDONE)) {
                        dx = tnv.ix - bp.ix;
                        dy = tnv.iy - bp.iy;
                        dz = tnv.iz - bp.iz;
                        square = dx * dx + dy * dy + dz * dz;
                        if (square < vpDistance[index]) {
                            boundPoint.set(tnv.ix, tnv.iy, tnv.iz, bp);

                            vpDistance[index] = square;
                            if (!(vpBits[index] & ISBOUND)) {
                                vpBits[index] |= ISBOUND;
                                outarray.push({
                                    ix : tnv.ix,
                                    iy : tnv.iy,
                                    iz : tnv.iz
                                });
                            }
                        }
                    }
                }
            }
        }

        // console.log("part3", positout);
        return outarray;
    };

    this.marchingcubeinit = function(stype) {
        for ( var i = 0, lim = vpBits.length; i < lim; i++) {
            if (stype == 1) {// vdw
                vpBits[i] &= ~ISBOUND;
            } else if (stype == 4) { // ses
                vpBits[i] &= ~ISDONE;
                if (vpBits[i] & ISBOUND)
                    vpBits[i] |= ISDONE;
                vpBits[i] &= ~ISBOUND;
            } else if (stype == 2) {// after vdw
                if ((vpBits[i] & ISBOUND) && (vpBits[i] & ISDONE))
                    vpBits[i] &= ~ISBOUND;
                else if ((vpBits[i] & ISBOUND) && !(vpBits[i] & ISDONE))
                    vpBits[i] |= ISDONE;
            } else if (stype == 3) { // sas
                vpBits[i] &= ~ISBOUND;
            }
        }
    };

    // this code allows me to empirically prune the marching cubes code tables
    // to more efficiently handle discrete data
    var counter = function() {
        var data = Array(256);
        for ( var i = 0; i < 256; i++)
            data[i] = [];

        this.incrementUsed = function(i, j) {
            if (typeof data[i][j] === 'undefined')
                data[i][j] = {
                    used : 0,
                    unused : 0
                };
            data[i][j].used++;
        };

        this.incrementUnused = function(i, j) {
            if (typeof data[i][j] === 'undefined')
                data[i][j] = {
                    used : 0,
                    unused : 0
                };
            data[i][j].unused++;

        };

        var redoTable = function(triTable) {
            var str = "[";
            for ( var i = 0; i < triTable.length; i++) {
                var code = 0;
                var table = triTable[i];
                for ( var j = 0; j < table.length; j++) {
                    code |= (1 << (table[j]));
                }
                str += "0x" + code.toString(16) + ", ";
            }
            str += "]";
            console.log(str);
        };

        this.print = function() {

            var table = MarchingCube.triTable;
            var str;
            var newtable = [];
            for ( var i = 0; i < table.length; i++) {
                var newarr = [];
                for ( var j = 0; j < table[i].length; j += 3) {
                    var k = j / 3;
                    if (typeof data[i][k] === 'undefined' || !data[i][k].unused) {
                        newarr.push(table[i][j]);
                        newarr.push(table[i][j + 1]);
                        newarr.push(table[i][j + 2]);
                    }
                    if (typeof data[i][k] === 'undefined')
                        console.log("undef " + i + "," + k);
                }
                newtable.push(newarr);
            }
            console.log(JSON.stringify(newtable));
            redoTable(newtable);
        };
    };
    
    this.marchingcube = function(stype) {
        this.marchingcubeinit(stype);
        verts = []; faces = [];   
        MarchingCube.march(vpBits, verts, faces, {
            smooth : 1,
            nX : pLength,
            nY : pWidth,
            nZ : pHeight        
        });      

        var pWH = pWidth*pHeight;
        for (var i = 0, vlen = verts.length; i < vlen; i++) {
            verts[i]['atomid'] = vpAtomID[verts[i].x * pWH + pHeight *
                    verts[i].y + verts[i].z];
        }  

        MarchingCube.laplacianSmooth(1, verts, faces);

    };


};


//generate a mesh from atoms
function generateMesh(atoms)
{
	var atomlist = []; //which atoms of atoms are selected for surface generation (all for us)
	$.each(atoms, function(i, a) { a.serial=i; atomlist[i]=i;}); //yeah, this is messed up 
	var time = new Date();
	
	var extent = $3Dmol.getExtent(atoms);
	var w = extent[1][0] - extent[0][0];
	var h = extent[1][1] - extent[0][1];
	var d = extent[1][2] - extent[0][2];
	var vol =  w * h * d;

	var ps = new ProteinSurface();
	ps.initparm(extent, true, vol);
	ps.fillvoxels(atoms, atomlist);
	ps.buildboundary();
	//compute solvent excluded
	ps.fastdistancemap();
	ps.boundingatom(false);
	ps.fillvoxelswaals(atoms, atomlist);	
	ps.marchingcube($3Dmol.SurfaceType.SES);

	var VandF = ps.getFacesAndVertices(atomlist);	
	
	//create geometry from vertices and faces
	var geo = new $3Dmol.Geometry(true);

	var geoGroup = geo.updateGeoGroup(0);

	var vertexArray = geoGroup.vertexArray;
	// reconstruct vertices and faces
	var v = VandF['vertices'];
	var offset;
	var i, il;
	for (i = 0, il = v.length; i < il; i++) {
		offset = geoGroup.vertices * 3;
		vertexArray[offset] = v[i].x;
		vertexArray[offset + 1] = v[i].y;
		vertexArray[offset + 2] = v[i].z;
		geoGroup.vertices++;
	}

	var faces = VandF['faces'];
	geoGroup.faceidx = faces.length;
	geo.initTypedArrays();

	// set colors for vertices
	var colors = [];
	for (i = 0, il = atoms.length; i < il; i++) {
		var atom = atoms[i];
		if (atom) {
			if (typeof (atom.surfaceColor) != "undefined") {
				colors[i] = atom.surfaceColor;
			} else if (atom.color) // map from atom
				colors[i] = $3Dmol.CC.color(atom.color);
		}
	}

	var verts = geoGroup.vertexArray;
	var colorArray = geoGroup.colorArray;
	var normalArray = geoGroup.normalArray;
	var vA, vB, vC, norm;

	// Setup colors, faces, and normals
	for (i = 0, il = faces.length; i < il; i += 3) {

		// var a = faces[i].a, b = faces[i].b, c = faces[i].c;
		var a = faces[i], b = faces[i + 1], c = faces[i + 2];
		var A = v[a]['atomid'];
		var B = v[b]['atomid'];
		var C = v[c]['atomid'];

		var offsetA = a * 3, offsetB = b * 3, offsetC = c * 3;

		colorArray[offsetA] = colors[A].r;
		colorArray[offsetA + 1] = colors[A].g;
		colorArray[offsetA + 2] = colors[A].b;
		colorArray[offsetB] = colors[B].r;
		colorArray[offsetB + 1] = colors[B].g;
		colorArray[offsetB + 2] = colors[B].b;
		colorArray[offsetC] = colors[C].r;
		colorArray[offsetC + 1] = colors[C].g;
		colorArray[offsetC + 2] = colors[C].b;

		// setup Normals

		vA = new $3Dmol.Vector3(verts[offsetA], verts[offsetA + 1],
				verts[offsetA + 2]);
		vB = new $3Dmol.Vector3(verts[offsetB], verts[offsetB + 1],
				verts[offsetB + 2]);
		vC = new $3Dmol.Vector3(verts[offsetC], verts[offsetC + 1],
				verts[offsetC + 2]);

		vC.subVectors(vC, vB);
		vA.subVectors(vA, vB);
		vC.cross(vA);

		// face normal
		norm = vC;
		norm.normalize();

		normalArray[offsetA] += norm.x;
		normalArray[offsetB] += norm.x;
		normalArray[offsetC] += norm.x;
		normalArray[offsetA + 1] += norm.y;
		normalArray[offsetB + 1] += norm.y;
		normalArray[offsetC + 1] += norm.y;
		normalArray[offsetA + 2] += norm.z;
		normalArray[offsetB + 2] += norm.z;
		normalArray[offsetC + 2] += norm.z;

	}
	geoGroup.faceArray = new Uint16Array(faces);
	var mat = new $3Dmol.MeshLambertMaterial();
	mat.vertexColors = true;
	
	var mesh = new $3Dmol.Mesh(geo, mat);
	mesh.doubleSided = true;
	
	var time2 = new Date();
	var t = time2-time;
	console.log("Surface mesh generation: "+t+"ms");
	$('#timeresult').html(t+"ms");
	return mesh;

}
