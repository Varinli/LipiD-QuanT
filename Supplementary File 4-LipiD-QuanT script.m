function lipidQuant = lipidQuant(filename)

%lipidQuant identifies and counts lipid droplets in cells
%
%   lipidQuant takes the filename of an image I as its input, and returns a
%   binary image BW of the same size as I, with 1's where the function
%   finds the droplets mask I and 0's elsewhere.
%
%   
%   Example
%   -------
%   Find lipid droplts in images D21-2-2.tif
%
%       I = lipidQuant('C:\D21-2-2.tif') 
%


%Image normalisation
close all
A=imread(filename);
B=double(A(:,:,1));
Bmean=mean(B(:));
bxstd=std2(B(:));
B=B-Bmean;B=(B*11.52)/bxstd;
B=B+101.69;
A=B;

% Laplacian edge detector.
temp=1-(edge(Gauss2D((double((A))),0.5),'zerocross',0));%


% Massaging the skeleton to avoid loosing droplets
temp=1-bwmorph(1-temp,'bridge');
temp(:,1)=0;temp(:,end)=0;temp(1,:)=0;temp(end,:)=0;
temp1=bwlabeln(temp,4);
temp1=imclose(1-(temp1>0),strel('square',1));
temp1=bwlabeln(1-temp1,4);
temp1(1:3,:)=max(temp1(:))+1;temp1((end-3):end,:)=max(temp1(:))+1;temp1(:,(end-3:end))=max(temp1(:))+1;temp1(:,1:3)=max(temp1(:))+1;

% Produce a mask temp1copy containing only droplet boundaries
temp1props=regionprops(temp1,del2(Gauss2D(medfilt2((double(A))),1)),'Area','PixelIdxList','MeanIntensity','Eccentricity');
temp1copy=zeros(size(temp1));
for k=1:size(temp1props,1)
    if temp1props(k).MeanIntensity>0.10                 
        temp1copy(temp1props(k).PixelIdxList)=1;
    end
    if temp1props(k).Area<200                           
        temp1copy(temp1props(k).PixelIdxList)=0;
    end
    
    if temp1props(k).Area<200 && temp1props(k).MeanIntensity>0.15 && temp1props(k).Area>20      
        temp1copy(temp1props(k).PixelIdxList)=1;
    end
end

% create mask cutitthere to cut large droplet pairs (dumbbells) into two parts
temp1copy1=1-bwareaopen(temp1copy,1000);
temp23=-bwdist(1-temp1copy1);
temp23=imimposemin(temp23,((MyLocalMax1(imhmax(-temp23,3)))>0).*temp1copy1);
cutitthere=((watershed(temp23)==0)).*temp1copy1;

%create a mask for the background temp1copyx
temp1props=regionprops(bwlabeln(1-temp1copy),'Area','PixelIdxList');
temp1copyx=zeros(size(temp1));
for k=1:size(temp1props,1)
    if temp1props(k).Area>10000
        temp1copyx(temp1props(k).PixelIdxList)=1;
    end
end

% Fish out small vesicles that would be lost otherwise
ta1=imopen(temp1copyx,strel('disk',1));
ta1bwn=bwlabeln(ta1);
tempa1props=regionprops(ta1bwn,'Area','PixelIdxList');
temp1copyxa3=zeros(size(temp1));

for k=1:size(tempa1props,1)
    if (size(tempa1props(k).PixelIdxList,1)<40)
        temp1copyxa3(tempa1props(k).PixelIdxList)=1;
    end
end
temp1copy(1:3,:)=1;temp1copy((end-3):end,:)=1;temp1copy(:,(end-3:end))=1;temp1copy(:,1:3)=1;

% create a mask for cell positions (no vesicle should be created outside mask)
mask=rangefilt(Gauss2D((double(A)),1))>1.8;
mask=imclose(mask,strel('disk',3));
temxxss=bwareaopen((1-mask),3000);
mask=bwareaopen((1-temxxss),3000);
mask=imerode(mask,strel('disk',5));
temp1copyxa3=temp1copyxa3.*mask; 
temp1copy=imopen(1-temp1copy,strel('disk',2));

% Only keep potential vesicles and fill any hole in vesicles
temp1=bwlabeln(temp1copy,4);
temp1props=regionprops(temp1,'All');

temp1copy=zeros(size(temp1));
for k=1:size(temp1props,1)
    if temp1props(k).Area<5000
        temp1copy(temp1props(k).PixelIdxList)=1;
    end
end
temp1copy=temp1copy.*mask;
temp1copy=imfill(temp1copy,'holes');

%Identify all large bits with negative curvature to eliminate void regions
%inside cells

temp1props=regionprops(bwlabeln(temp1copy,4),'All');
temp1copynoq=zeros(size(temp1));
for k=1:size(temp1props,1)
    try
        if temp1props(k).Area>200
            temp1copynoq(temp1props(k).PixelIdxList)=1;
        else continue
        end
        if  temp1props(k).Area>200&&((temp1props(k).Area/temp1props(k).ConvexArea)>0.86)
            temp1copynoq(temp1props(k).PixelIdxList)=0;
            continue
        end       
        % Measure curvature 
        patch=zeros(200);
        coos=round(temp1props(k).PixelList-repmat(temp1props(k).Centroid-[100,100],[size(temp1props(k).PixelList,1),1]));
        patch(sub2ind([size(patch,1),size(patch,2)],coos(:,1),coos(:,2)))=1;
        patch=Gauss2D(patch,1);patch=imresize(patch,3);
        patch=patch>0.5;
        [listCONTOUR,listNORMALS] = TRACE_MooreNeighbourhood(patch);%TRACE_MooreNeighbourhood was contributed by Adam H. Aitkenhead 
        listCONTOUR4=listCONTOUR;
        for i=4:(size(listCONTOUR,1)-4)
            listCONTOUR4(i,:)=mean(listCONTOUR(i-3:i+3,:));
        end       
        Vertices=listCONTOUR4(1:4:end,:)*10;
        Lines=[(1:size(Vertices,1))' (2:size(Vertices,1)+1)']; Lines(end,2)=1;
        km=LineCurvature2D(Vertices,Lines);
        km=km*10000;              
        if size(find(km>0),1)>(size(km,1)/2.4)
            temp1copynoq(temp1props(k).PixelIdxList)=1;
        else
            temp1copynoq(temp1props(k).PixelIdxList)=0;
        end       
    catch
        continue
    end
end

% Split dumbbell droplets if any
temp1copy=temp1copy.*(1-(cutitthere));

% Final quality control on droplets

temp1props=regionprops(bwlabeln(temp1copy),'All');
temp1copynoq=zeros(size(temp1));
for k=1:size(temp1props,1)
    if temp1props(k).Area<5000
        temp1copynoq(temp1props(k).PixelIdxList)=1;
    end
    if temp1props(k).Area>100 && temp1props(k).Eccentricity>0.90
        temp1copynoq(temp1props(k).PixelIdxList)=0;
    end
    if  temp1props(k).Area>200&&((temp1props(k).Area/temp1props(k).ConvexArea)<0.60)
        temp1copynoq(temp1props(k).PixelIdxList)=0;
    end
end


% Produce binary result image in same directory as input image
imwrite(temp1copynoq,[filename(1:end-4),'out.tiff'],'tiff');

% Export a text file containing the area of every droplet in the image
lipidQuant=(temp1copynoq+temp1copyxa3)>0;
temp=bwlabeln(temp1copynoq);
tempo=regionprops(temp);
fid = fopen([filename(1:end-4),'out.txt'],'w');
a='%12.2f ';
fprintf(fid,a,([tempo.Area]));
fclose(fid);



