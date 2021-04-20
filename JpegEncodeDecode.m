clear all

% Q = csvread('luminance.csv');
Q = [16 11 10 16 24 40 51 61; %luminance
    12 12 14 19 26 58 60 55;
    14 13 16 24 40 57 69 56;
    14 17 22 29 51 87 80 62;
    18 22 37 56 68 109 103 77;
    24 35 55 64 81 104 113 92;
    49 64 78 87 103 121 120 101;
    72 92 95 98 112 100 103 99];

% Fungsi 1 - Membaca gambar dan mengubah RGB menjadi Greyscale
imdata = imread('fruits.jpg');
[rows, columns,color] = size(imdata);

G = rgb2gray(imdata);
Size=size(G);
num_mem_bytes = prod(Size)

clear imdata;

% Fungsi 2 - Mengubah gambar menjadi matrix
columnsnew = columns;
rowsnew = rows;

if mod(rows,8) ~= 0
    rowsnew = rows + mod(rows,8);
end

if mod(columns,8) ~= 0
    columnsnew = columns + mod(columns,8);
end

A = zeros(rowsnew,columnsnew);


for i =1:1:columns
    for j =1:1:rows
        A(j,i) = G(j,i);
    end
end

clear columnsnew rowsnew;

% Fungsi 3 - Membagi matriks menjadi 8x8
for i=1:8:columns - 8
    for j=1:8:rows - 8
        newrows = (j+7)/8;
        newcolumns = (i+7)/8;
        B{newrows,newcolumns} = A(j:j+7, i:i+7);
    end
end

clear newrows newcolumns;

% Fungsi 4 - Membangun matrix DCT
DCT = ones(8,8);

for g = 0:7
    for j = 0:7
        DCT(j+1,g+1) = cos(((2*g + 1)*j*pi)/16);
    end
end

% ENCODING

for i = 1:149
    for j = 1:78
        B{j,i} = DCT * B{j,i} * transpose(DCT) ;
        
        % Fungsi 5 - Mengkuantisasi matrix dengan luminance
        L = round(B{j,i}/(Q));
        
        
        % Fungsi 6 - Melakukan Run Length Encoding
        X=reshape(L',[1,64]);
        M=rle(X);
        
        
        % Fungsi 7 - Melakukan Huffman Encoding
        symbols = unique(M(:));
        counts = hist(M(:), symbols);
        p = double(counts) ./ sum(counts);
        sig = M;
        dict{j,i} = huffmandict(symbols,p);
        H{j,i} = huffmanenco(sig,dict{j,i});
    end
end



% Print Encoded File
HTable = cell2table(H);
dictTable = cell2table(dict);
writetable(HTable, 'EncodeResult.csv');
writetable(dictTable, 'EncodeDictionary.csv');

% clearvars -except Q H dict

% DECODING

for i = 1:149
    for j = 1:78
        % Fungsi 8 - Melakukan Huffman Decoding
        dsig = huffmandeco(H{j,i}', dict{j,i});
        
        % Fungsi 9 - Inverse Run Length Encoding
        J = irle(dsig');
        J = reshape(J,[8,8]);
        J = J';
        
        % Fungsi 10 - Dekuantisasi matrix
        D = round(J*Q);
        
        % Fungsi 11 - Melakukan Operasi Inverse DCT
        hasil{j,i} = DCT\D/transpose(DCT);
    end
end

% Fungsi 12 - Menampilkan gambar
hasil = cell2mat(hasil);
hasil = uint8(hasil);
imshow(hasil);
title('After Encoding and Decoding');
imwrite(hasil, 'fruitsdecoded.jpg');
M = Size(1);
N = Size(2);
newHasil = zeros(size(G));
newHasil(1:size(hasil,1),1:size(hasil,2)) = hasil;
hasil2 = uint8(newHasil);
MSE = sum(sum((G-hasil2).^2))/(M*N);
PSNR = 10*log10(256*256/MSE);

% Fungsi Lanjutan
% Fungsi 6 - Melakukan Run Length Encoding
function Output=rle(Input)
Len=length(Input);
v=1;
k=1;
o=1;
while v<2*Len
    comp=1;
    for v=v:Len
        if v==Len
            break
        end
        if Input(v)==Input(v+1)
            comp=comp+1;
        else
            break
        end
    end
    Output(k)=comp;
    Output(k+1)=Input(v);
    if v==Len && Input(v-1)==Input(v)
        break
    end
    o=o+1;
    k=k+2;
    v=v+1;
    if v==Len
        if mod(Len,2)==0
            Output(k)=1;
            Output(k+1)=Input(v);
        else
            Output(k)=1;
            Output(k+1)=Input(v);
        end
        break
    end
end
end

% Fungsi 9 - Melakukan Decode dari RLE
function Output=irle(Input)
Len=length(Input);
newlength = 0;
v = 0;
b = 0;
s = 0;
x = [];
for v=1:2:Len-1
    for b = 1:1:Input(v)
        x = [x Input(v+1)];
    end
end
Output=x;
end