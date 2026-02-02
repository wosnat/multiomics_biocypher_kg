function gsea2(expfile,rankfile,varargin)
% Performs a variant of GSEA as described in Subramanian et.al. PNAS (early
% ed -- DOI 10.1073.  rvals are the 'correlation' values for each gene, and
% they must be in sort order.  svals are +1 or -1 depending on whether the
% gene is in the reference set S.  pexp is the exponent for the rvals for
% computation of the ES value.  Gene permutation only is used here, so
% there is no need for the actual data (there are not enough conditions in
% the data I am using to make condition permutation workable).  No multiple
% hypothesis testing is supported.
%
% This program, gsea2, is a version of gsea1 that is set up to process
% Prochlorococcus data given me by Andy Tolonen on 10/17/2005.  The
% procedure now reads a file and does only one form of gsea test, rather
% than several forms as done by gsea1.
% 
% Required parameters:
%  expfile = name of input expression data file in format produced by Andy
%  rankfile = rankings of NtcA sites for genes as given by Su.et.al.(2005)
%    NAR.  File consists of three columns: First column is gene name,
%    second is rank as in Su.et.al. supplementary data, third column is 
%    relative rank for all genes with orthologs, 4th column is relative
%    rank for genes without orthologs.  Subsequent columns are ignored.
%    Following format of Su et.al. data, it is possible for multiple genes
%    to be registered on a single file line (if they are in a predicted
%    operon).  In this case, the genes are given the same ranks. The option 
%    to recognize ortholog vs. non-ortholog rankings is supported here 
%    because Su et.al. rankings are heavily influenced by presence of orthologs.
%
% Keyword parameters: given as keyword, value pairs     
%   'titl': character string used in titles of figures
%   'numperms': number of permutations used to compute enrichment score p-value
%      Default = 5000 (1000 used by per Subramanian et.al. (2005)).
%   'pexp': exponent used to compute ES score.  Default = 1 as per
%      Subramanian et.al. 2005.
%   'OverallRank': highest rank in overall Su et.al. (2005) NtcA scores to
%      be considered part of 'S' set tested for enrichment in expression
%      data.  Default = 25.  'OverallRank' is not used if either
%      'OrthologRank' or 'NonOrthologRank' is provided.
%   'OrthologRank': highest relative rank of orthologs as per Su et.al NtcA
%      scores to be considered part of 'S' set tested for enrichment in
%      expression data.  Default = 0 (not used).
%   'NonOrthologRank': highest relative rank of non-orthologs as per Su
%      et.al. NtcA scores to be considered part of 'S' set tested for
%      enrichment in expression data.  Default = 0 (not used).
%   'outpref': prefix for output file name.  Default = expfile + ".out"
%     
%
% Created: 10/17/2005-10/18/2005 by John Aach
% Copyright (c) 2005 by John Aach and the President and Fellows of Harvard
% University

global parameters rvals svals genelist generanks
parameters=strvcat('titl','pexp','numperms','outpref', ...
    'OverallRank', 'OrthologRank','NonOrthologRank');
rvals=[];
svals=[];
generanks=[];

for i=1:size(parameters,1)
    currparm=deblank(parameters(i,:));
    evalstr=sprintf('global %s %s_d;',currparm,currparm);
    eval(evalstr);
    evalstr=sprintf('%s=[];',currparm);
    eval(evalstr);
end;
titl_d=[];
pexp_d=1;
numperms_d=5000;
outpref_d=[];
OverallRank_d=25;
OthologRank_d=0;
NonOrthologRank_d=0;

if nargin<1
    error('Input expression file not provided');
end;
if nargin<2
    error('Input NtcA rank file not provided');
end;

parmstr=AnalyzeParameters(expfile,varargin{:});
outfile=[outpref,'.txt'];

ReadExpFile(expfile);
GetGeneRanks(rankfile);

% get distribution of es scores based on randomization
for i=1:numperms
    es_perm=compute_ES(rvals,svals(randperm(length(svals))),pexp);
    es_dist(i)=es_perm;
end;

[es_actual,leadedge_actual,es_scores_actual]=compute_ES(rvals,svals,pexp);
if es_actual>max(es_dist(:))
    pval=1/numperms;
else
    pval=sum(es_dist(:)>=es_actual)/numperms;
end;

% Chebyschev bound
es_dist_std=std(es_dist(:));
es_dist_mean=mean(es_dist(:));
es_actual_Z = abs(es_actual-es_dist_mean)/es_dist_std;
pval_cheby=1/(es_actual_Z.^2);

% print output file
outfid=fopen(outfile,'w');
fprintf(outfid,'PARAMETERS FOR THIS EXECUTION OF gsea2.m:\n%s\nEND OF PARAMETERS\n\n',parmstr);
fprintf(outfid,'Max ES score:\t%g\n',es_actual);
fprintf(outfid,'P value:\t%g\n',pval);
fprintf(outfid,'Null ES mean:\t%g\n',es_dist_mean);
fprintf(outfid,'Null ES std:\t%g\n',es_dist_std);
fprintf(outfid,'Max ES score Z score:\t%g\n',es_actual_Z);
fprintf(outfid,'P value Chebyshev:\t%g\n',pval_cheby);
fprintf(outfid,'Leading Edge Size:\t%d',leadedge_actual);
fprintf(outfid,'\n\nGene\tES_score\tInLeadingEdge?\tInSet?\tOverallRank\tOrthologRank\tNonOrthologRank\n');
for i=1:size(genelist,1)
    fprintf(outfid,'%s\t%g\t%d\t%d\t%d\t%d\t%d\n',deblank(genelist(i,:)),es_scores_actual(i), ...
        i<=leadedge_actual,svals(i)==1,generanks(i,1),generanks(i,2),generanks(i,3));
end;
fclose(outfid);

% plot figures
figure;
plot(1:length(es_scores_actual),es_scores_actual(:));
yy=ylim;
xx=xlim;
hold on;
plot([leadedge_actual leadedge_actual],[yy(1) yy(2)],'r:');
PrintData(xx,yy,es_actual,pval,pval_cheby,leadedge_actual);
title([titl,': score plot']);
xlabel('gene number');
ylabel('running enrichment score');

figure;
[n,xout]=hist(es_dist(:));
bar(xout,n);
hold on;
yy=ylim;
plot([es_actual es_actual],[yy(1) yy(2)],'r:');
xx=xlim;
PrintData(xx,yy,es_actual,pval,pval_cheby);
title([titl,': null distribution']);
xlabel('gene number');
ylabel('running enrichment score');

return;

function PrintData(xx,yy,es_value,pval,pval_cheby,leadedge)
  h_es=text(mean(xx(:)),mean(yy(:)),sprintf('ES=%g',es_value),'FontSize',8);
  extent_es=get(h_es,'Extent');
  width_es=extent_es(3);
  height_es=extent_es(4);
  h_P=text(mean(xx(:)),mean(yy(:)),sprintf('P=%g',pval),'FontSize',8);
  extent_P=get(h_P,'Extent');
  width_P=extent_P(3);
  height_P=extent_P(4);
  max_width=max(width_es,width_P);
  max_height=max(height_es,height_P);
  h_Pc=text(mean(xx(:)),mean(yy(:)),sprintf('P_c=%g',pval_cheby),'FontSize',8);
  extent_Pc=get(h_P,'Extent');
  width_Pc=extent_P(3);
  height_Pc=extent_P(4);
  max_width=max(max_width,width_Pc);
  max_height=max(max_height,height_Pc);
  if nargin==6
      h_le=text(mean(xx(:)),mean(yy(:)),sprintf('le=%d',leadedge),'FontSize',8);
      extent_le=get(h_le,'Extent');
      width_le=extent_le(3);
      height_le=extent_le(4);
      max_width=max(max_width,width_le);
      max_height=max(max_height,height_le);
  end;
  set(h_es,'Position',[xx(2)-max_width yy(2)-max_height]);
  set(h_P,'Position',[xx(2)-max_width yy(2)-2*max_height]);
  set(h_Pc,'Position',[xx(2)-max_width yy(2)-3*max_height]);
  if nargin==6
      set(h_le,'Position',[xx(2)-max_width yy(2)-4*max_height]);
  end;
  return;

function [es,leadedge_ix,es_scores_all]=compute_ES(rvals,svals,pexp);
  nargout_num=nargout;
  N_hit=0;
  norm_hit_score=zeros(length(rvals),1);
  norm_hit=zeros(length(rvals),1);
  for i=1:length(rvals);
      if svals(i)==1
          N_hit=N_hit+1;
          norm_hit(i)=1;
          if pexp~=0
              if rvals(i)~=0
                  norm_hit_score(i)=abs(rvals(i))^pexp;
              else
                  norm_hit_score(i)=0;
              end;
          else
              norm_hit_score(i)=1;
          end;
      end;
  end;
  norm_hit_score=cumsum(norm_hit_score);
  norm_hit_score=norm_hit_score/norm_hit_score(end);
  N_miss=length(rvals)-N_hit;
  norm_miss_score=cumsum(1-norm_hit)/N_miss;
  
  es_score=zeros(length(rvals),1);
  for i=1:length(rvals)
      es_score(i)=norm_hit_score(i)-norm_miss_score(i);      
      if i==1
          max_es_score_ix=i;
      elseif es_score(i)>es_score(max_es_score_ix)
          max_es_score_ix=i;
      end;
  end;
  if nargout_num>=1
      es=es_score(max_es_score_ix);
  end;
  if nargout_num>=2
      leadedge_ix=max_es_score_ix;
  end;
  if nargout_num>=3
      es_scores_all=es_score;
  end;
  
  return;
      
function parmstr=AnalyzeParameters(expfile,varargin);
  global parameters 
  
  for i=1:size(parameters,1)
      currparm=deblank(parameters(i,:));
      evalstr=sprintf('global %s %s_d;',currparm,currparm);
      eval(evalstr);
  end;
  
  parmprovided=zeros(size(parameters,1),1);
  i=1;
  while(i<=length(varargin));
      currkeyword=varargin{i};
      match_ix=0;
      matchparm=[];
      for j=1:size(parameters,1);
          currparm=deblank(parameters(j,:));
          if ~isempty(regexpi(currparm,['^',currkeyword]))
              if match_ix==0
                  match_ix=j;
                  matchparm=currparm;
              else
                  error(sprintf('Ambiguous keyword %s matches multiple parameters.',currkeyword));
              end;
          end;
      end;
      if match_ix==0
          error(sprintf('Unrecognized keyword %s',currkeyword));
      end;
      parmprovided(match_ix)=1;
      if strcmpi(matchparm,'titl') | ...
              strcmpi(matchparm,'outpref') | ...
              strcmpi(matchparm,'pexp') | ...
              strcmpi(matchparm,'numperms') | ...
              strcmpi(matchparm,'OverallRank') | ...
              strcmpi(matchparm,'OrthologRank') | ...
              strcmpi(matchparm,'NonOrthologRank')
          if i==length(varargin)
              error(sprintf('Value for keyword parameter %s not provided.', matchparm));
          end;
          currval=varargin{i+1};
          if strcmpi(matchparm,'titl')
              if ~ischar(currval)
                  error(sprintf('Invalid value provided for parameter %s',matchparm));
              end;
              titl=currval;
              i=i+1;
          elseif strcmpi(matchparm,'outpref')
              if ~ischar(currval)
                  error(sprintf('Invalid value provided for parameter %s',matchparm));
              end;
              outpref=currval;
              i=i+1;
          elseif strcmpi(matchparm,'pexp')
              if ~isnumeric(currval) 
                  error(sprintf('Invalid value provided for parameter %s',matchparm));
              end;
              pexp=currval;
              i=i+1;
          elseif strcmpi(matchparm,'numperms')
              if ~isnumeric(currval) | currval<=0 | floor(currval) ~= currval 
                  error(sprintf('Invalid value provided for parameter %s',matchparm));
              end;
              numperms=currval;
              i=i+1;
          elseif strcmpi(matchparm,'OverallRank')
              if ~isnumeric(currval) | currval<=0 | floor(currval) ~= currval 
                  error(sprintf('Invalid value provided for parameter %s',matchparm));
              end;
              OverallRank=currval;
              i=i+1;
          elseif strcmpi(matchparm,'OrthologRank')
              if ~isnumeric(currval) | currval<=0 | floor(currval) ~= currval 
                  error(sprintf('Invalid value provided for parameter %s',matchparm));
              end;
              OrthologRank=currval;
              i=i+1;
          elseif strcmpi(matchparm,'NonOrthologRank')
              if ~isnumeric(currval) | currval<=0 | floor(currval) ~= currval 
                  error(sprintf('Invalid value provided for parameter %s',matchparm));
              end;
              NonOrthologRank=currval;
              i=i+1;
          else 
              error(sprintf('Unprocessable parameter %s',matchparm));
          end;
      end;
      i=i+1;
  end;
          
  if (OrthologRank>0 | NonOrthologRank>0) & length(OverallRank)>0 
      error('Inconsistent specification: OverallRank incompatible with OrthologRank and NonOrthologRank');
  end;
  for i=1:length(parmprovided)
      if ~parmprovided(i)
          currparm=deblank(parameters(i,:));
          evalstr=sprintf('%s=%s_d;',currparm,currparm);
          eval(evalstr);
      end;
  end;
  if (OrthologRank>0 | NonOrthologRank>0)
      OverallRank=0;
  end;
  
  if length(outpref)==0
      if length(expfile)>=4 & strcmpi(expfile(end-4:end),'.txt')
          outpref=[expfile(1:end-4),'.out'];
      else
          outpref=[expfile,'.out'];
      end;
  end;
  
  parmstr=[];
  if length(titl)>0
      parmstr=sprintf('%stitl=''%s''',parmstr,titl);
  end;
  
  if length(outpref)>0
      if length(parmstr)>0
          parmstr=sprintf('%s\n',parmstr);
      end;
      parmstr=sprintf('%soutpref=''%s''',parmstr,outpref);
  end;
  
  if length(parmstr)>0
      parmstr=sprintf('%s\n',parmstr);
  end;
  parmstr=sprintf('%snumperms=%d',parmstr,numperms);
  
  if length(parmstr)>0
      parmstr=sprintf('%s\n',parmstr);
  end;
  parmstr=sprintf('%spexp=%g',parmstr,pexp);

  if OverallRank>0
      if length(parmstr)>0
          parmstr=sprintf('%s\n',parmstr);
      end;
      parmstr=sprintf('%sOverallRank=%d',parmstr,OverallRank);
  end;
     
  if OrthologRank>0
      if length(parmstr)>0
          parmstr=sprintf('%s\n',parmstr);
      end;
      parmstr=sprintf('%sOrthologRank=%d',parmstr,OrthologRank);
  end;

  if NonOrthologRank>0
      if length(parmstr)>0
          parmstr=sprintf('%s\n',parmstr);
      end;
      parmstr=sprintf('%sNonOrthologRank=%d',parmstr,NonOrthologRank);
  end;
  return;
  
function [gene,data]=ParseLine(fileline);
% parses file line that has a row header (gene) + data
 [gene,rem]=strtok(fileline,9);  % 9 = tab
 i=0;
 data={};
 while(length(rem)>0)
     [value,rem]=strtok(rem,9);
     i=i+1;
     if length(value)==0
         data{i}=[];
         continue;
     end;
     [numval,count,errmsg,nextix]=sscanf(value,'%g');  % convoluted logic, but isnumeric doesn't seem to work...
     if length(errmsg)==0
         data{i}=numval;
     else
         data{i}=value;
     end;
 end;
 return;
  
function ReadExpFile(expfile);
  global rvals genelist
  rvals=[];
  genelist=[];
  
  infid=fopen(expfile);
  if infid==0
      error(sprintf('Could not open expression data file %s',expfile));
  end;
  
  dataline=fgetl(infid);  % skip title line
  signed1mpval=[];
  
  while(~feof(infid))
      dataline=fgetl(infid);
      [gene,data]=ParseLine(dataline);
      genelist=strvcat(genelist,gene);
      i=size(genelist,1);
      if sign(data{1})~=0
          signed1mpval(i,1:2)=[i sign(data{1})*(1-exp(data{3}))];
      else
          signed1mpval(i,1:2)=[i (1-2.^data{3})];
      end;
  end;
 fclose(infid);
 signed1mpval=sortrows(signed1mpval,2);
 signed1mpval=signed1mpval(end:-1:1,:);
 rvals=signed1mpval(:,2);
 genelist=genelist(signed1mpval(:,1),:);
 return;

function GetGeneRanks(rankfile);
  global svals genelist generanks

  generanks=zeros(size(genelist,1),3);
  
  infid=fopen(rankfile);
  if infid==0
      error(sprintf('Could not open gene rank file %s',rankfile));
  end;
  
  dataline=fgetl(infid);  % skip title line
  svals=-ones(size(genelist,1),1);
  
  while(~feof(infid))
      dataline=fgetl(infid);
      [gene,data]=ParseLine(dataline);
      
      % parse 'gene', as this may contain multiple whitespace delimited gene
      % names
      [gene1,rem]=strtok(gene);
      RegisterGeneInSvals(gene1,data{:})
      while(length(rem)>0)
          [gene1,rem]=strtok(rem);
          RegisterGeneInSvals(gene1,data{:});
      end;
  end;
  fclose(infid);
  return;
      
function RegisterGeneInSvals(gene,varargin)         
  global svals genelist OrthologRank NonOrthologRank OverallRank generanks
  data=varargin;
  i=1;
  foundgene=0;
  while(i<=size(genelist,1) & ~foundgene)
      currgene=deblank(genelist(i,:));
      if strcmpi(currgene,gene)
          foundgene=1;
      else
          i=i+1;
      end;
  end;
  if ~foundgene
      return;
  end;
  
  generanks(i,1)=data{1};  % save data
  generanks(i,2)=data{2};
  generanks(i,3)=data{3};
  
  if OverallRank>0 & data{1}>0 & data{1}<=OverallRank
      svals(i)=1;
  else
      if OrthologRank>0 & data{2}>0 & data{2}<=OrthologRank
          svals(i)=1;
      end;
      if NonOrthologRank>0 & data{3}>0 & data{3}<=NonOrthologRank
          svals(i)=1;
      end;
  end;
  return;