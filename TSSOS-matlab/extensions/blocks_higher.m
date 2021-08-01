function [data] = blocks_higher(data)
%BLOCKS_HIGHER apply the block closure operation again for a tssos psatz

fbasis=data.fbasis;
gbasis=data.gbasis;
% fsizes=data.fsizes;
% fsupp=data.fsupp;
% gsupp=data.gsupp;
ssupp=data.ssupp;
% coe=data.coe;
lt=data.lt;
fblocks = data.fblocks;
gblocks = data.gblocks;

fblocksize = data.fblocksize;
gblocksize = data.gblocksize;

m = length(data.gbasis);
[flb, n] = size(data.fbasis);

glb = cellfun(@(s) size(s, 1), data.gbasis);

%% extract the support of monomials inside the blocks
fcl=length(fblocks);
gcl = cellfun(@length, gblocks);
    fsupp=[];
    for i=1:fcl
        for j=1:fblocksize(i)
            for r=j:fblocksize(i)
                bi=fbasis(fblocks{i}(j),:)+fbasis(fblocks{i}(r),:);
                fsupp=[fsupp;bi];
            end
        end
    end
    gsupp=cell(1,m);
%     supp1=fsupp;
    for k=1:m
        gsupp{k}=[];
        for i=1:gcl(k)
            for j=1:gblocksize{k}(i)
                for r=j:gblocksize{k}(i)
                    for s=1:lt(k+1)
                        bi=ssupp{k+1}(s,:)+gbasis{k}(gblocks{k}{i}(j),:)+gbasis{k}(gblocks{k}{i}(r),:);
                        gsupp{k}=[gsupp{k};bi];
                    end
                end
            end
        end
%         supp1=[supp1;gsupp{k}];
    end
    
    
%% perform a new block closure operation
ofsupp=odd_supp(n,fsupp);
ofsupp=unique(ofsupp,'rows');
lofsize=size(ofsupp);
lfo=lofsize(1);
flb=length(fbasis);
fedges=[];
for i = 1:flb
    for j = i:flb
        bi=fbasis(i,:)+fbasis(j,:);
        if length(bi(~mod(bi,2)))==n||bfind(ofsupp,lfo,bi,n)~=0
           fedges=[fedges;[i j]];
        end
    end
end
fG=graph(fedges(:,1),fedges(:,2));
fblocks=conncomp(fG,'OutputForm','cell');
fcl=length(fblocks);
fblocksize=zeros(1,fcl);
for i=1:fcl
    fblocksize(i)=length(fblocks{i});
end
% nfsizes=[unique(fblocksize);hist(fblocksize,unique(fblocksize))];
[gc,grps] = groupcounts(fblocksize);
nfsizes = flip([grps, gc], 1)';
fsize_old = data.fblocksize;
s1=size(fsize_old);
s2=size(nfsizes);             
if s1(2)~=s2(2)||~isempty(setdiff(nfsizes,fsize_old,'rows'))
   data.fsizes=nfsizes;
   disp('fblocksizes=');
   disp(nfsizes);
else
   opt=0;
   data=0;
   status=0;
   disp('No higher blocking hierarchy'); 
   return
end
disp('gblocksizes=');
glb=zeros(1,m);
gblocks=cell(1,m); 
gcl=zeros(1,m);
gblocksize=cell(1,m);
for k=1:m
    glb(k)=length(gbasis{k});
    ggsupp=[fsupp;gsupp{k}];
    ogsupp=odd_supp(n,ggsupp);
    ogsupp=unique(ogsupp,'rows');
    logsize=size(ogsupp);
    lgo=logsize(1);
    gedges=[];
       for i = 1:glb(k)
           for j = i:glb(k)
               r=1;
               while r<=lt(k+1)
                     bi=ssupp{k+1}(r,:)+gbasis{k}(i,:)+gbasis{k}(j,:);
                     if length(bi(~mod(bi,2)))==n||bfind(ogsupp,lgo,bi,n)~=0
                        break
                     else
                        r=r+1;
                     end
               end
               if r<=lt(k+1)
                  gedges=[gedges;[i j]];
               end
           end
       end
       gG=graph(gedges(:,1),gedges(:,2));
       gblocks{k}=conncomp(gG,'OutputForm','cell');
       gcl(k)=length(gblocks{k});
       gblocksize{k}=zeros(1,gcl(k));
       for i=1:gcl(k)
           gblocksize{k}(i)=length(gblocks{k}{i});
       end
%     gsizes=[unique(gblocksize{k});hist(gblocksize{k},unique(gblocksize{k}))];
    [gc,grps] = groupcounts(gblocksize{k});
    gsizes = flip([grps, gc], 1)';
    disp(gsizes);
end
% supp1=unique(supp1,'rows');



%% store the data
data.fblocks=fblocks;
data.gblocks=gblocks;

data.fblocksize=fblocksize;
data.gblocksize=gblocksize;

end

