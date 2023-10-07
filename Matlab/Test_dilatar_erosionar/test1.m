BW = imread('1600_1200_1.jpg')
se = strel('disk', 3, 8);
Bw1_dilatado = imdilate(BW,se);
Bw1_erosionado = imerode(BW,se);
BW_D_E = imerode(Bw1_dilatado,se);
BW_E_D = imerode(Bw1_erosionado,se);

%%C = imfuse(A,B) crea una imagen compuesta a partir de dos imágenes, 
% A y B. Si A y B son de diferente tamaño, imfuse rellena con ceros las 
% dimensiones más pequeñas para que las dos imágenes tengan el mismo tamaño. 
% La salida, C, es una matriz numérica 
% que contiene una versión fusionada de las imágenes A y B.

%figure, imshow(Bw1_dilatado), title('Dilatado');
%figure, imshow(Bw1_erosionado), title('Erosionado');
subplot(3,2,1), imshow(Bw1_dilatado), title('Dilatado') ;
subplot(3,2,2), imshow(Bw1_erosionado), title('Erosionado') ;
subplot(3,2,3), imshow(BW_D_E), title('Dilatado - Erosionado') ;
subplot(3,2,4), imshow(BW_E_D), title('Erosionado - Dilatado') ;
subplot(3,2,5), imshow(BW), title('Normal') ;
%figure, imshow(BW)