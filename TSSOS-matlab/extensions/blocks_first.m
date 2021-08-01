function [data] = blocks_first(data)
%BLOCKS_FIRST Find PSD blocks for the first closure operation in a tssos
%psatz. Takes in data from get_support_monom.m

%Adds the following fields to data:
%   fsizes: fsizes
%   fsizes: fsupp
%   fsizes: gsupp


%% load in the previously gathered data
fbasis=data.fbasis;
gbasis=data.gbasis;
% fsizes=data.fsizes;
% fsupp=data.fsupp;
% gsupp=data.gsupp;
ssupp=data.ssupp;
% coe=data.coe;
lt=data.lt;

m = length(data.gbasis);
[flb, n] = size(data.fbasis);

glb = cellfun(@(s) size(s, 1), data.gbasis);
%% copied from blockpop_cons_first.m 

%find the support preserved among all polynomials
supp = vertcat(ssupp{:});

supp=unique(supp,'rows', 'sorted');



%obtain the blocks
osupp=odd_supp(n,supp);
losize=size(osupp);
lo=losize(1);
fedges=[];
for i = 1:flb
    for j = i:flb
        bi=fbasis(i,:)+fbasis(j,:);
        %if (all exponents are even) or (sum of exponents are in
        %odd-support basis)
        if length(bi(~mod(bi,2)))==n||bfind(osupp,lo,bi,n)~=0
           fedges=[fedges;[i j]];
        end
    end
end
fG=graph(fedges(:,1),fedges(:,2));
fblocks=conncomp(fG,'OutputForm','cell');
fcl=length(fblocks);
if fcl==1
   opt=0;
   data=0;
   status=0;
   disp('Unblockable');
   return
end
fblocksize=zeros(fcl,1);
for i=1:fcl
    fblocksize(i)=length(fblocks{i});
end
% fsizes=[unique(fblocksize);hist(fblocksize,unique(fblocksize))];
[gc,grps] = groupcounts(fblocksize);
fsizes = flip([grps, gc], 1)';
disp('fblocksizes=');
disp(fsizes);
disp('gblocksizes=');
gblocks=cell(1,m); 
gcl=zeros(1,m);
gblocksize=cell(1,m);
for k=1:m
    gedges=[];  
    for i = 1:glb(k)
        for j = i:glb(k)
            r=1;
            while r<=lt(k+1)
                  bi=ssupp{k+1}(r,:)+gbasis{k}(i,:)+gbasis{k}(j,:);
                  if length(bi(~mod(bi,2)))==n||bfind(osupp,lo,bi,n)~=0
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
    gblocksize{k}=zeros(gcl(k),1);
    for i=1:gcl(k)
        gblocksize{k}(i)=length(gblocks{k}{i});
    end
    [gc,grps] = groupcounts(gblocksize{k});
    gsizes = flip([grps, gc], 1)';
%     gsizes=[unique(gblocksize{k});hist(gblocksize{k},unique(gblocksize{k}))];
    disp(gsizes);
end

%% find supports in the block decomposition



%% store the data
data.fblocks=fblocks;
data.gblocks=gblocks;

data.fblocksize=fblocksize;
data.gblocksize=gblocksize;
% data.fsupp = fsupp;
% data.gsupp = gsupp;
% data.supp1 = supp1;

%when do fsupp and gsupp get defined?
%there should be a way to do this without calling blockcpop and solving the
%sdp. This was fixed in Julia

%[fsupp,gsupp,opt,status]=blockcpop(n,m,ssupp,coe,lt,fbasis,gbasis,fblocks,fcl,fblocksize,gblocks,gcl,gblocksize,solver);
end

