function [DayNum, TideHeight] = ReadWATideData(FileName, RefYear)
%This function reads "txt" historical tide data downloaded from the WA
%Department of Transport website:
%https://www.transport.wa.gov.au/imarine/download-tide-wave-data.asp
%
%Usage: 
%[DayNum, TideHeight] = ReadWATideData(FileName, RefYear);
%FileName is a string containing the name of the file to be read (including path and extension)
%RefYear - A DayNum of 0 corresponds to 00:00:00 on the first of January of
%this year (numeric)
%
%DayNum is a column vector of Matlab date numbers (see help on datenum function)
%TideHeight is a corresponding column vector of tide heights in metres
%
%Alec Duncan, 5/4/2023

DateFmt = 'yyyymmdd.HHMM';
BadCm = -9999;  %Bad data value

Fid = fopen(FileName, 'rt');
if Fid <= 0
    error(['Couldn''t open file: ' FileName]);
end

disp(['Reading: ' FileName]);

%Read the file header and display
GotHeader = false;
while ~GotHeader
    InStr = fgetl(Fid);
    if strcmp(InStr(1), '#')
        disp(InStr);
    else
        GotHeader = true;
        %Skip this line as it is the header for the data values
    end
end


AllDatStr = fscanf(Fid, '%s');  %Read the rest of the data as one giant string
fclose(Fid);

AllDat = strsplit(AllDatStr, {',', '/'}); %Break it up into individual strings for height and date/time

NAll = floor(length(AllDat)/2) * 2; 
TideHeight = (str2double(AllDat(1:2:NAll)).')/100;  %Odd elements ar eheights
DayNum = datenum(AllDat(2:2:NAll), DateFmt) - datenum([RefYear, 1, 1, 0, 0, 0]);  %Even elements are date/times

%Delete any bad data values
IsBad = abs(TideHeight - BadCm/100) < eps;
TideHeight(IsBad) = [];
DayNum(IsBad) = [];



