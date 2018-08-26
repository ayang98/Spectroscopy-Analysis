function RamanAnalysis
%ExcelTransfer
PeakFinder(300,400)%1511, 1292 vx-770 ACN 
%correction(1511)
%plot_spec(150,1890) %150- 1890 is full wavelength range

end

function correction(peak) %IMPORTANT- use this on each group of 32 spectra to get 32 different points at certain peaks corrected
[filename, pathname] = uigetfile('.spc', 'MultiSelect', 'on');
option = input('Do you have one SPC file (enter 1) or multiple SPC files (enter 2)? : ');
[wave,specs,fid] = Spectroscopy(filename,pathname,option);
keep = 1:1600;
MS = mscorrect(specs);
plot(specs(:,1890-peak), 'r*')
title('Uncorrected')
xlabel('Spectra Number')
ylabel('Intensity (counts)')
figure();
plot(MS(:,1890-peak), 'r*')
title('MS Corrected VX-770 SDD at 1511 cm^-1')
xlabel('Time (min)')
ylabel('Intensity (counts)')
end


function plot_spec(range1,range2)
[filename, pathname] = uigetfile('.spc', 'MultiSelect', 'on');
option = input('Do you have one SPC file (enter 1) or multiple SPC files (enter 2)? : ');
[wave,specs,fid] = Spectroscopy(filename,pathname,option);
ind = wave > range1 & wave < range2;
w1 = wave(ind);
s1 = specs(:,ind);
PlotAllSpectra(w1,s1,filename)
xlabel('Raman shift (cm^{-1})')
ylabel('Intensity (counts)')
legend 
end

function PeakFinder(range1, range2) %IMPORTANT- can only use for one spectra at a time (out of possible 32 for each experiment)
[filename, pathname] = uigetfile('.spc', 'MultiSelect', 'on');
%option = input('Do you have one SPC file (enter 1) or multiple SPC files (enter 2)? : ');
[wave,specs,fid] = Spectroscopy(filename,pathname,1); %only used currently for one spectra at a time1
PlotAllSpectra(wave,specs,filename)
figure();
[w1, s1] = trimspectra(wave,specs,range1, range2);
%disp(w1)
hold on
[pp] = cubicspline(w1,s1);
[z, d1] = intercepts(w1,pp);
%disp(z)
interpolate(z,w1, d1, pp, s1);
%for i=1:length(z)
    %disp(w1(z(i)));
%end
end

function interpolate(z, w1, d1, pp, s1)
f = @(xx) -ppval(pp,xx); %multipy by negative one so that the mins of this new function are the maxs of the old function
%since matlab only has fminbnd

% now loop through each value found in z to find the minimum over the
% interval around each zero.
x = zeros(size(z));
for i=1:length(x)
    if z(i)+1<=length(w1)
        lowerbound = w1(z(i)+1);
    else
        lowerbound = w1(z(i));
    end
    if z(i)-1 >= 1
        upperbound = w1(z(i)-1);
    else
        upperbound = w1(z(i));
    end
    x(i) = fminbnd(f,lowerbound,upperbound);
end
figure
plot(w1, ppval(d1,w1))
hold on
plot(x, 0,'ro')

xlabel('Raman shift (cm^{-1})')
legend('1^{st} derivative', 'zeros')

figure(2)
plot(x, ppval(pp,x), 'ro','markerfacecolor','g')

% print some offset labels

%disp(length(x))
for i=1:length(x)
     out = sprintf('%1.1f cm^{-1}',x(i));
     disp(strcat([num2str(i),':',' ',out])); %square brackets for spaces
 
end

end

function [z, d1] = intercepts(w1,pp)
d1 = fnder(pp, 1); % first derivative
s = ppval(d1, w1);
s(s >= 0) = 1;
s(s < 0) = 0;
% where diff = -1 indicates a transition from positive slope to
% negative slope, which is a maximum. These are approximate locations
% of the zeros
%d=diff(s)
z = find(diff(s) == 1); %IMPORTANT- need to set equal to one because you have flipped the direction of x axis!!!! 
end


function [pp] = cubicspline(w1, s1)
pp = csaps(w1, s1,.99);
hold all
plot(w1, ppval(pp,w1))
set(gca,'xdir','reverse');
xlabel('Raman shift (cm^{-1})')
ylabel('Intensity (counts)')
end


function [w1,s1] = trimspectra(wave,specs,range1, range2)
ind = wave > range1 & wave < range2;
w1 = wave(ind);
s1 = specs(ind);
figure(2)
plot(w1, s1,'r*')
set(gca,'xdir','reverse');
xlabel('Raman shift (cm^{-1})')
ylabel('Intensity (counts)')

end



function [max_shift] = maxSpec(wave,specs)
M = max(specs);
for i= 1:length(specs)
    if specs(i)==M
        max_shift = wave(i);
        break
    end
end
end


function ExcelTransfer %This function asks for path of folder storing Raman data as well as a name
%of the newly created excel file as input and returns an excel file with
%all the Raman data
tic
prompt1 = 'Pathname of Raman data folder: ';
in1 = input(prompt1,'s');
prompt2 = 'Desired Excel file name: ';
in2 = input(prompt2,'s');
prompt3 = 'Pathname of MATLAB folder: ';
in3 = input(prompt3,'s');
A = dir(in1); %A is a struct
cell = {A.name}; %get name field as cell array
for i=3:length(cell)-1 %trim cell array (get rid of first 2 and last weird things)
    names(i)=cell(i);
end
names = names((3:end));
%disp(names)
filename2 = strcat(in2,'.xlsx'); %change name of newly created Excel file here 
for i=1:length(names)%need to loop over all folders, then all spectra in each folder (2 loops)
    s1= strcat(in1,'\');
    s2= names{i}; %do this because cell array {}
    total = strcat(s1,s2);
    
    %disp(total);
    C = dir(total);
    cell = {C.name};
    %disp(cell)
    for a=3:length(cell)
        names2(a)=cell(a);
    end
    %disp(names2)
    names2 = names2((3:end)); %important because first two are erroneous 
   
    filename = {names2{1},names2{16},names2{length(names2)}};
    pathname = strcat(total,'\');
    [wave,specs,fid] = Spectroscopy(filename,pathname,2);
    keep = 1:1600;
    %PlotAllSpectra(wave(keep),specs(:,keep)) %check to make sure doing the right thing (71
    %windows, 3 plots each window)
    
    sheet = i;
    for j = 1:3
        xlRange = strcat('A',num2str(3*j-2));
        xlswrite(filename2,[wave(keep); specs(j,keep)],sheet,xlRange);
    end
    %{
    xlRange = 'A1';
    xlswrite(filename2,[wave(keep); specs(1,keep)],sheet,xlRange);
    xlRange = 'A4';
    xlswrite(filename2,[wave(keep); specs(2,keep)],sheet,xlRange);
    xlRange = 'A7';
    xlswrite(filename2,[wave(keep); specs(3,keep)],sheet,xlRange);
    %}
end
e = actxserver('Excel.Application'); % # open Activex server
ewb = e.Workbooks.Open(strcat(in3,'\',in2)); % # open file (enter full path!)
for i=1:length(names)
    ewb.Worksheets.Item(i).Name = names{i};
end
ewb.Save % # save to the same file
ewb.Close(false)
e.Quit
toc
end




function secPlot(specs)
keep = 1:1600;
MS = mscorrect(specs(:,keep));
T = linspace(0,32,32);
plot(T, MS(:,280), 'r*')
[data,p] = secder(specs(:,keep), 2, 15);
figure();
plot(data')
end



function PlotAllSpectra(wave, specs, filename)
S = size(specs);
first_dim = S(1);
if first_dim ==1
    %plot(specs);
    %plot(fliplr(specs(1:length(wave))));
    %plot(specs(length(wave):length(specs))); %%not sure why but specs vector not trimmed
    %to length of waves vector for just one graph
    i=1;
    %figure();
    plot(wave,specs(i,:));
    set(gca,'xdir','reverse'); %reverse the horizontal orientation of axis so goes from highest to lowest wavelength 
    
else
    i=1;
    figure();
    hold on 
    while(i<first_dim+1)
        plot(wave,specs(i,:));
        set(gca,'xdir','reverse'); %reverse the horizontal orientation of axis so goes from highest to lowest wavelength 
        i=i+1;
    end
end
xlabel('Raman shift (cm^{-1})')
ylabel('Intensity (counts)')
if first_dim == 1 %creates the title- for once spec file, title is filename which is simply a string
    new = '';
    number = '';
    for j = length(filename)-5:length(filename)-4
        number = strcat(number, filename(j));
    end
    %disp(number)
    for i=1:length(filename)-8 %minus 9 to get rid of '.spc' and long number
        if filename(i)~='_'
            new = strcat(new,filename(i));
        end
    end  
    disp(new)
    if ismember(new(length(new)),['N','L']) == 0
        new2 = '';
        for i = 1:length(new)-1
            new2 = strcat(new2,new(i));
        end
        title(strcat(new2,number))
    else
        title(strcat(new,number))
    end
else %for multiple spec files, title is the first item in filename, which is a cell array
    first_file = filename{1};
    new = '';
    for i=1:length(first_file)-8 %minus 8 to get rid of number and '.spc'
        if first_file(i)~='_'
            new = strcat(new,first_file(i));
        end
    end  
    %disp(new)
    if ismember(new(length(new)),['N','L']) == 0
        new2 = '';
        for i = 1:length(new)-1
            new2 = strcat(new2,new(i));
        end
        title(new2)
    else
    title(new)
    end
end
disp('Plot(s) generated');
end

function [wave,specs,fid] = Spectroscopy(filename,pathname,option)    %RamanRxn2

%[filename, pathname] = uigetfile('.spc', 'MultiSelect', 'on');
disp(pathname)
%option = input('Do you have one SPC file (enter 1) or multiple SPC files (enter 2)? : ');


wave = 150:1890;
%% For Multiple SPC files:
if option == 2;
    
    [~,nc] = size(filename);
 
% Import specs    
    for i = 1:nc
%         fid(i,:) = cell2mat(filename(i));
        a = fopen([pathname cell2mat(filename(i))]);
        b = fread(a,inf,'float');
        specs(i,:) = b(wave);
        fclose('all');
    end
% Determine Time and tElapsed for individual SPC files...
   ... Updated 20Aug2010
fid = [filename, pathname]';
% [nr,nc] = size(fid);
% 
%     for i = 1:(nc-1);
%         file = filename(i);
%         a = [pathname file];
%         b = cell2mat(a);
%         c = fopen([b]);
%         d = fread(c, '*char')';
%            
%  %       hour = str2num(d(14008:14009));
%  %   if isempty(hour)== 0;
%  %       min = str2num(d(14011:14012));
%  %       sec = str2num(d(14014:14015)); 
%  %   else
%  %        hour = str2num(d(14009:14010));
%  %        min = str2num(d(14012:14013));
%  %        sec = str2num(d(14015:14016));
%  %   end
% 
%    
%  %     Time(i) = ((60*hour) + (min) + (sec/60));
%  %       tElapsed(i) = Time(i) - Time(1);
%    
%         fclose all; 
%    end
%  
%% If one SPC file:
elseif option == 1

    a = fopen([pathname filename]);
    b = fread(a,inf,'float');
    c = b(137:end);
    nospecs = floor(length(c)/3335);
    for i = 1:nospecs
        %specs(i,:) = c(1:3327);
        specs(i,:) = b(wave); %specs needs to be same length as waves
        c = c(3335:end);
    end
end
fid = [filename, pathname]';
%% Wavelength
% wave = 150:1890;
wave = fliplr(wave);
%fid = [pathname filename];

%fid = [pathname filename];

%if option == 2;
%    [tElapsed, Time] = TimeElapsedSPC(fid);
%elseif option == 1;
end

function outspec = mscorrect(tnspec)
[nrows,ncols] = size(tnspec);
outspec = zeros(nrows,ncols);
%ones(ncols,1)
%tnspec'
%mean(tnspec)' %take the mean of tnspec first, then transpose that IMPORTANT- mean(A) takes mean of COLUMNS if A is matrix!!!\
if nrows ==1
    X = [ones(ncols,1),tnspec'];
    disp(X)
else
    X = [ones(ncols,1),mean(tnspec)'];
end
for i=1:nrows
    coeff = regress(tnspec(i,:)',X);
    outspec(i,:) = (tnspec(i,:)-coeff(1))/coeff(2);
end

end

function [data,p] = secder(inspec, order, tol_sg)

%% inspec is the vector of spectra
%% order is the order, e.g. 1 for first derivative, 2 for second
%% derivative
%% This function calculates the second derivative using a polynomial fit.
% tol_sg = savitzky-golay filter width (i.e. 5,7,9,11)

a = 200;
b = -3.611;
tol = a*tol_sg^b;

%tol = 0.001;
[nrows,ncols] = size(inspec);
%x_ax = [1:ncols];
data = zeros(nrows,ncols);
if nrows <= 100;   
    %[pp,p] = csaps([1:ncols],inspec,0.0000001);
    %[pp,p] = csaps(x_ax,inspec);
    [pp,p] = csaps(1:ncols,inspec,tol);
    [ppder] = fnder(pp,order);
    data = fnval(ppder,1:ncols);
elseif nrows > 100;    
    h = waitbar(0,'Initializing Waitbar');
    for i = 1:nrows;
         temp = i/nrows;
         if mod(i,50) == 0;
             waitbar(temp,h,sprintf('%2.1f%% done..',100*temp));
         elseif mod(i,50) ~= 0;
         end
        %[pp,p] = csaps([1:ncols],inspec(i,:));
        [pp,p] = csaps(1:ncols,inspec(i,:),tol);
        [ppder] = fnder(pp,order);
        data(i,:) = fnval(ppder,1:ncols);        
    end
    close(h);
end

end