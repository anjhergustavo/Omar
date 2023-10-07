%Rejilla con dimensiones 1600x1200
x = linspace(-0.01, 0.01, 1600); %distancias en metros, 1600 pix in X
T = 0.8E-3 %Periodo de la rejilla de difraccion
y = cos(2*pi*(1/T)*x);
y = y >= 0;
y = repmat(y,1200,1); %Las dimensiones se hacen con linspace y repmat
imshow(y); 