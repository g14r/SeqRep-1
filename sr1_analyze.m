function[varargout]=sr1_analyze(what,varargin)
%% [varargout]=sr1_analyze(what,varargin)
% SequenceRepetition experiment (sr1) analysis of behavioral data
% usage|example calls:
%                       sr1_analyze('all_subj');                            %pre-analysis: run the subject routine for all_subj
%                       sr1_analyze('all_subj',{'s01'});                    %pre-analysis: run the subject routine for selected subjects
%                       sr1_analyze('all_data');                            %pre-analysis: create .mat file containing data from subjects in subj
%                       D=sr1_analyze('repetitionEffect');                  %group analysis: repetition effect
%                       D=sr1_analyze('repetitionEffect',4);                %group analysis: repetition effect (collapsing repNum>=4)
%                       D=sr1_analyze('repetitionEffect','indiv_subj');     %single subject analysis: repetition effect separately in each subject
%                       D=sr1_analyze('repetitionEffect','indiv_subj',4);   %single subject analysis: repetition effect separately in each subject (collapsing repNum>=4)
%                       D=sr1_analyze('repetitionEffect','s01');            %single subject analysis: repetition effect in subject s01
%                       D=sr1_analyze('repetitionEffect','s01',4);          %single subject analysis: repetition effect in subject s01 (collapsing repNum>=4)
%                       sr1_analyze('pointsAnalysis');                      %proportion of points/errors (effect of thresholds) per block number,subject and group
%
%                       D=sr1_analyze('repEff_ET');
%                       D=sr1_analyze('repEff_RT');
%                       D=sr1_analyze('repDiff');
%                       D=sr1_analyze('repCorr');
%                       D=sr1_analyze('repSff');
% --
% gariani@uwo.ca - 2017.07.17

%% paths
pathToData='/Users/gariani/Documents/data/SequenceRepetition/sr1';
pathToAnalyze='/Users/gariani/Documents/data/SequenceRepetition/sr1/analyze';

%% globals
subj={'s01','s02','s03','s04','s05','s06','s07','s08','s09','s10',...
      's11',      's13','s14','s15','s16','s17','s18','s19','s20',...
      's21','s22','s23','s24','s25','s26','s27'}; %'s04' is the experimenter, 's12' was excluded due to not following instructions
ns = numel(subj);
subvec = zeros(1,ns);
for i = 1:ns; subvec(1,i) = str2double(subj{i}(2:3)); end

% plot defaults
fs = 20; %default font size for all figures
lw = 4; %4; %default line width for all figures
ms = 10; %10; %12; %default marker size for all figures
rs = 0;%0:4; %[0,1,2,3,4]; %default repetition number subset (full range = 0:4; 0 = full set)
maxRep = 0; %3; % default ceiling level for repetition number (0 = no ceiling)

% colors
cbs_red = [213 94 0]/255;
cbs_blue = [0 114 178]/255;
cbs_yellow = [240 228 66]/255;
cbs_pink = [204 121 167]/255;
cbs_green = [0 158 115]/255;
blue = [49,130,189]/255;
lightblue = [158,202,225]/255;
red = [222,45,38]/255;
lightred = [252,146,114]/255;
green = [49,163,84]/255;
lightgreen = [161,217,155]/255;
orange = [253,141,60]/255;
yellow = [254,196,79]/255;
lightyellow = [255,237,160]/255;
purple = [117,107,177]/255;
lightpurple = [188,189,220]/255;
darkgray = [50,50,50]/255;
gray = [150,150,150]/255;
lightgray = [200,200,200]/255;
silver = [240,240,240]/255;
black = [0,0,0]/255;
white = [255,255,255]/255;

% styles
style.reset;
style.custom({blue,lightblue,red,lightred,orange,yellow,lightyellow,purple,lightpurple,darkgray,gray,lightgray,green,lightgreen,black,silver,white,...
    cbs_red,cbs_yellow,cbs_blue,cbs_green,cbs_pink});
isrepsty = style.custom({lightgray, darkgray}, 'markersize',ms, 'linewidth',lw, 'errorbars','shade');
lightgraysty = style.custom({lightgray}, 'markertype','none', 'linewidth',1);
graysty = style.custom({lightgray}, 'markersize',ms, 'linewidth',lw, 'errorbars','shade');
darkgraysty = style.custom({darkgray}, 'markersize',ms, 'linewidth',lw, 'errorbars','shade');
blacksty = style.custom({black}, 'markertype','none', 'linewidth',lw, 'linestyle','--','errorbars','plusminus', 'errorwidth',lw, 'errorcap',0);%, 'linestyle','none');
%sffsty = style.custom({lightgray, gray, darkgray}, 'markersize',ms, 'linewidth',lw);
%boxplotsty = style.custom({darkgray}, 'markertype','none', 'linewidth',lw);
isrepstybox = style.custom({gray, black}, 'markersize',lw, 'linewidth',lw);
daysty = style.custom({gray,darkgray}, 'markersize',ms, 'linewidth',lw);

cr={[1 0 0],[.9 .9 .9],[.7 .7 .7],[.5 .5 .5],[.3 .3 .3]}; %default line color per repNum
cd={[.2 .6 .5],[.6 .5 .8]}; %default line color per day

% legends
isrepleg = {'Switch', 'Repetition'};
dayleg = {'Day 1','Day 2'}; %default legend per day
subl = {'s01','s01','s02','s02','s03','s03','s05','s05','s06','s06','s07','s07','s08','s08','s09','s09','s10','s10',...
    's11','s11','s13','s13','s14','s14','s15','s15','s16','s16','s17','s17','s18','s18','s19','s19','s20','s20'};
seql = [1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8];

%% types of analysis
switch (what)
    case 'all_subj' % pre-analysis: run the subject routine for all_subj
        if nargin>1; subj=varargin{1}; end
        for s=1:numel(subj)
            sr1_subj(subj{s},0); %run sr1_subj.m routine (without plot)
        end
        
    case 'all_data' % pre-analysis: create .mat file containing data from all subjects
        all_data=[];
        for s=1:length(subj)
            fprintf('\n%s\n\n',subj{s});
            D=load(fullfile(pathToData,sprintf('sr1_%s.mat',subj{s}))); %load data structure for each subject
            D.SN=ones(numel(D.TN),1)*s; %add information about subject number
            D.day=ones(numel(D.TN),1); D.day(D.BN>12)=2; %add information about day of testing
            D.repNum = D.repNum-1; % 0 = switch; 1 = first repetition (i.e., sequence repeated twice in a row)
            D.isRep = D.repNum>0;
            D = getrow(D, D.dummy == 0);
            D.ET = D.MT;
            
            %-------------------------------------------------------------------------------------------------------------------------------------
            % add IPI info
            D.IPI = diff([D.pressTime1, D.pressTime2, D.pressTime3, D.pressTime4], 1, 2);
            D.IPI_1 = D.IPI(:, 1); D.IPI_2 = D.IPI(:, 2); D.IPI_3 = D.IPI(:, 3);
            
            %-------------------------------------------------------------------------------------------------------------------------------------
            B = [];
            for b = 1:max(D.BN)
                T = getrow(D, D.BN==b);
                %-------------------------------------------------------------------------------------------------------------------------------------
                % detect Same First Finger (sff) in different sequences
                T.sff = zeros(numel(T.TN), 1);
                for t = 2:numel(T.sff)
                    if (T.isRep(t) == 1)
                        % repetition
                        T.sff(t) = 2;
                    elseif (T.isRep(t) == 0) && (T.press1(t, 1) == T.press1(t-1, 1))
                        % switch, same first finger (sff)
                        T.sff(t) = 1;
                    elseif (T.isRep(t) == 0)
                        % switch
                        T.sff(t) = 0;
                    else % not categorized
                        error('This transition does not fall into any of the defined categories!');
                    end
                end
                %-------------------------------------------------------------------------------------------------------------------------------------
                % add Same Finger Trasnistion info
                T.sft = zeros(numel(T.TN), 1);
                for t = 2:numel(T.sft)
                    % build transition matrix for consecutive sequences
                    seq1 = num2str(T.cueP(t));
                    tm1 = zeros(5);
                    for i = 2:(numel(seq1))
                        tm1(str2double(seq1(i)), str2double(seq1(i-1))) = 1;
                    end
                    seq2 = num2str(T.cueP(t-1));
                    tm2 = zeros(5);
                    for i = 2:(numel(seq2))
                        tm2(str2double(seq2(i)), str2double(seq2(i-1))) = 1;
                    end
                    if (T.isRep(t) == 1)
                        % repetition
                        T.sft(t) = 123;
                    elseif (T.isRep(t) == 0) && (T.press2(t, 1) == T.press2(t-1, 1)) && (T.press3(t, 1) == T.press3(t-1, 1)) && (T.press4(t, 1) == T.press4(t-1, 1))
                        % switch, same scond and third transitions
                        T.sft(t) = 23;
                    elseif (T.isRep(t) == 0) && (T.press1(t, 1) == T.press1(t-1, 1)) && (T.press2(t, 1) == T.press2(t-1, 1)) && (T.press3(t, 1) == T.press3(t-1, 1))
                        % switch, same first and second transitions
                        T.sft(t) = 12;
                    elseif (T.isRep(t) == 0) && (T.press3(t, 1) == T.press3(t-1, 1)) && (T.press4(t, 1) == T.press4(t-1, 1))
                        % switch, same third transition
                        T.sft(t) = 3;
                    elseif (T.isRep(t) == 0) && (T.press2(t, 1) == T.press2(t-1, 1)) && (T.press3(t, 1) == T.press3(t-1, 1))
                        % switch, same second transition
                        T.sft(t) = 2;
                    elseif (T.isRep(t) == 0) && (T.press1(t, 1) == T.press1(t-1, 1)) && (T.press2(t, 1) == T.press2(t-1, 1))
                        % switch, same first transition
                        T.sft(t) = 1;
                    elseif (T.isRep(t) == 0) && ~any(any(tm1==1 & tm1==tm2))
                        % switch, all transitions are different
                        T.sft(t) = 0;
                    elseif (T.isRep(t) == 0)
                        % switch, some transitions may be shared
                        T.sft(t) = -1;
                    else
                        % not categorized
                        error('This transition does not fall into any of the defined categories!');
                    end
                end
                B = addstruct(B, T);
            end
            
            all_data = addstruct(all_data, B); %append data structures from each subject
        end
        save(fullfile(pathToAnalyze,'sr1_all_data.mat'),'-struct', 'all_data'); %save all_data.mat file
        
    case 'repetitionEffect' %repetition effect: group and single subject analysis
        % single subject analysis
        if nargin>1 && ~isnumeric(varargin{1})
            if strcmp(varargin(1),'indiv_subj'); else; subj=varargin(1); end %either for each and every subject in the group, or just one in particular
            for s=1:length(subj)
                D=load(fullfile(pathToData,sprintf('sr1_%s.mat',subj{s}))); %load data for this subject
                D.SN=ones(numel(D.TN),1)*s; %add information about subject number
                D.day=ones(numel(D.TN),1); D.day(D.BN>12)=2; %add information about day of testing
                if nargin>2; D.repNum(D.repNum>varargin{2})=varargin{2}; end %potentially collapse data for repNum>varargin{2}
                % create summary table for RT and MT
                D1=tapply(D,{'day','repNum'},...
                    {D.RT(:),'mean','name','RT'},...
                    {D.MT(:),'mean','name','MT'},...
                    'subset',D.dummy==0 & D.isError==0);
                % build figure legend per repNum
                reps=sort(unique(D1.repNum)); cr=cr(1:numel(reps)); lr=cell(1,numel(reps));
                for ir=1:numel(reps)
                    if ir==numel(reps) && numel(reps)<5; lr{1,ir}=sprintf('Rep %d+',reps(ir)); else; lr{1,ir}=sprintf('Rep %d',reps(ir)); end
                end
                % RT
                figure('Name',sprintf('Repetition effect - %s',subj{s}));
                subplot(3,4,1); title('Reaction Time (RT)');
                lineplot(D1.repNum,D1.RT,'split',D1.day,'linewidth',lw,'linecolor',cd,'leg',dayleg,'leglocation','northeast');
                set(gca,'xticklabel',lr); xlabel('Repetition number'); ylabel('Mean RT (ms)'); set(gca,'fontsize',fs);
                % MT
                subplot(3,4,2); title('Movement Time (MT)');
                lineplot(D1.repNum,D1.MT,'split',D1.day,'linewidth',lw,'linecolor',cd,'leg',dayleg,'leglocation','northeast');
                set(gca,'xticklabel',lr); xlabel('Repetition number'); ylabel('Mean MT (ms)'); set(gca,'fontsize',fs);
                % MT diff
                D2=tapply(D,{'day'},...
                    {D.MT,'mean','name','MTfirst','subset',D.repNum==1},...
                    {D.MT,'mean','name','MTrest','subset',D.repNum>1},...
                    'subset',D.dummy==0 & D.isError==0);
                subplot(3,4,[3,4]); title('First vs Rest'); hold on;
                barplot(D2.day,D2.MTfirst-D2.MTrest);
                set(gca,'xticklabel',[1 2]); xlabel('Day'); ylabel('MT diff (ms)'); set(gca,'fontsize',fs); drawline(0,'dir','horz','linewidth',2);
                % MT learning
                D3=tapply(D,{'BN','repNum'},...
                    {D.MT(:),'mean','name','MT'},...
                    'subset',D.dummy==0 & D.isError==0);
                subplot(3,4,[5,6]); title('Learning: Movement time (MT)');
                lineplot(D3.BN,D3.MT,'split',D3.repNum,'errorbars','plusminus','errorcap',0,'linewidth',lw,'linecolor',cr,'leg',lr,'leglocation','northeast'); drawline(12.5,'dir','vert','linewidth',10,'color',[1 1 1]); drawline(12.5,'dir','vert','linewidth',lw);
                xlabel('Block number'); ylabel('Mean MT (ms)'); set(gca,'fontsize',fs);
                D4=tapply(D,{'BN'},...
                    {D.MT,'mean','name','MTfirst','subset',D.repNum==1},...
                    {D.MT,'mean','name','MTrest','subset',D.repNum>1},...
                    'subset',D.dummy==0 & D.isError==0);
                subplot(3,4,[9,10]); title('MT First - MT Rest');
                lineplot(D4.BN,D4.MTfirst-D4.MTrest,'errorbars','shade','shadecolor',[.5 .5 .5],'linewidth',lw,'linecolor',[.5 .5 .5],'leg',lr,'leglocation','northeast'); drawline(0,'dir','horz','linewidth',lw,'linestyle','--'); drawline(12.5,'dir','vert','linewidth',10,'color',[1 1 1]); drawline(12.5,'dir','vert','linewidth',lw);
                xlabel('Block number'); ylabel('MT diff (ms)'); set(gca,'fontsize',fs);
                % error learning
                D5=tapply(D,{'BN','repNum'},...
                    {D.isError,'sum','name','sumErr'},...
                    {D.TN,'numel','name','nTrials'},...
                    'subset',D.dummy==0);
                D5.ER=((D5.sumErr./D5.nTrials)*100);
                subplot(3,4,[7,8]); title('Learning: Error rate (ER)');
                lineplot(D5.BN,D5.ER,'split',D5.repNum,'errorbars','plusminus','errorcap',0,'linewidth',lw,'linecolor',cr,'leg',lr,'leglocation','northeast'); drawline(15,'dir','horz','linewidth',lw,'linestyle','--');
                xlabel('Block number'); ylabel('Mean ER (%)'); set(gca,'fontsize',fs); set(gca,'ylim',[0 100]); drawline(12.5,'dir','vert','linewidth',10,'color',[1 1 1]); drawline(12.5,'dir','vert','linewidth',lw);
                D6=tapply(D5,{'BN'},...
                    {D5.ER,'mean','name','ERfirst','subset',D5.repNum==1},...
                    {D5.ER,'mean','name','ERrest','subset',D5.repNum>1});
                subplot(3,4,[11,12]); title('ER First - ER Rest');
                lineplot(D6.BN,D6.ERfirst-D6.ERrest,'errorbars','shade','shadecolor',[.5 .5 .5],'linewidth',lw,'linecolor',[.5 .5 .5],'leg',lr,'leglocation','northeast'); drawline(0,'dir','horz','linewidth',lw,'linestyle','--');
                xlabel('Block number'); ylabel('ER diff (%)'); set(gca,'fontsize',fs); drawline(12.5,'dir','vert','linewidth',10,'color',[1 1 1]); drawline(12.5,'dir','vert','linewidth',lw);
                
                % out
                D.D1=D1; D.D2=D2; D.D3=D3; D.D4=D4; D.D5=D5; D.D6=D6; %incorporate the sub-structures as fields of main structure
                varargout={D}; %return main structure for this subject
            end;
        else
            % group analysis
            D=load(fullfile(pathToAnalyze,'sr1_all_data.mat')); %load group data
            if nargin>1; D.repNum(D.repNum>varargin{1})=varargin{1}; end %potentially collapse data for repNum>varargin{1}
            % create summary table for RT and MT
            D1=tapply(D,{'SN','day','repNum'},...
                {D.RT(:),'mean','name','RT'},...
                {D.MT(:),'mean','name','MT'},...
                'subset',D.dummy==0 & D.isError==0);
            % build figure legend per repNum
            reps=sort(unique(D1.repNum)); cr=cr(1:numel(reps)); lr=cell(1,numel(reps));
            for ir=1:numel(reps)
                if ir==numel(reps) && numel(reps)<5; lr{1,ir}=sprintf('Rep %d+',reps(ir)); else; lr{1,ir}=sprintf('Rep %d',reps(ir)); end
            end
            % RT
            figure('Name',sprintf('Repetition effect - group (N=%d)',numel(subj)));
            subplot(3,4,1);
            lineplot(D1.repNum,D1.RT,'split',D1.day,'linewidth',lw,'linecolor',cd,'leg',dayleg,'leglocation','northeast');
            set(gca,'xticklabel',lr); xlabel('Repetition number'); ylabel('Mean RT (ms)'); set(gca,'fontsize',fs);
            [t3,p3]=ttest(D1.RT(D1.repNum==1&D1.day==1),D1.RT(D1.repNum==2&D1.day==1),2,'paired');
            [t5,p5]=ttest(D1.RT(D1.repNum==1&D1.day==2),D1.RT(D1.repNum==2&D1.day==2),2,'paired');
            title(sprintf('RT Rep1 vs Rep2:   tD1=%1.2f pD1=%1.2f   tD2=%1.2f pD2=%1.2f',t3,p3,t5,p5));
            % MT
            subplot(3,4,2);
            lineplot(D1.repNum,D1.MT,'split',D1.day,'linewidth',lw,'linecolor',cd,'leg',dayleg,'leglocation','northeast');
            set(gca,'xticklabel',lr); xlabel('Repetition number'); ylabel('Mean MT (ms)'); set(gca,'fontsize',fs);
            [t4,p4]=ttest(D1.MT(D1.repNum==1&D1.day==1),D1.MT(D1.repNum==2&D1.day==1),2,'paired');
            [t6,p6]=ttest(D1.MT(D1.repNum==1&D1.day==2),D1.MT(D1.repNum==2&D1.day==2),2,'paired');
            title(sprintf('MT Rep1 vs Rep2:   tD1=%1.2f pD1=%1.2f   tD2=%1.2f pD2=%1.2f',t4,p4,t6,p6));
            % MT diff
            D2=tapply(D,{'SN','day'},...
                {D.MT,'mean','name','MTfirst','subset',D.repNum==1},...
                {D.MT,'mean','name','MTrest','subset',D.repNum>1},...
                'subset',D.dummy==0 & D.isError==0);
            [t1,p1]=ttest(D2.MTfirst(D2.day==1)-D2.MTrest(D2.day==1),0,2,'onesample');
            [t2,p2]=ttest(D2.MTfirst(D2.day==2)-D2.MTrest(D2.day==2),0,2,'onesample');
            subplot(3,4,[3,4]); title(sprintf('Delta MT Rep1-RepOther vs 0:     tD1=%1.2f  pD1=%1.2f     tD2=%1.2f  pD2=%1.2f',t1,p1,t2,p2));
            myboxplot(D2.day,D2.MTfirst-D2.MTrest,'plotall',2,'xtickoff'); hold on; lineplot(D2.day,D2.MTfirst-D2.MTrest,'split',D2.SN,'linewidth',lw,'linecolor',[.5 .5 .5]);
            set(gca,'xticklabel',[1 2]); xlabel('Day'); ylabel('Mean Delta MT (ms)'); set(gca,'fontsize',fs); drawline(0,'dir','horz','linewidth',lw,'linestyle','--');
            % MT learning curves
            D3=tapply(D,{'SN','BN','repNum'},...
                {D.MT(:),'mean','name','MT'},...
                'subset',D.dummy==0 & D.isError==0);
            subplot(3,4,[5,6]); title('Learning: Movement time (MT)');
            lineplot(D3.BN,D3.MT,'split',D3.repNum,'errorbars','plusminus','errorcap',0,'linewidth',lw,'linecolor',cr,'leg',lr,'leglocation','northeast');
            drawline(12.5,'dir','vert','linewidth',28,'color',[1 1 1]); drawline(12.5,'dir','vert','linewidth',lw);
            xlabel('Block number'); ylabel('Mean MT (ms)'); set(gca,'fontsize',fs);
            D4=tapply(D,{'SN','BN'},...
                {D.MT,'mean','name','MTfirst','subset',D.repNum==1},...
                {D.MT,'mean','name','MTrest','subset',D.repNum>1},...
                'subset',D.dummy==0 & D.isError==0);
            subplot(3,4,[9,10]); title('Delta MT Rep1-RepOther');
            lineplot(D4.BN,D4.MTfirst-D4.MTrest,'errorbars','shade','shadecolor',[.5 .5 .5],'linewidth',lw,'linecolor',[.5 .5 .5],'leg',lr,'leglocation','northeast'); drawline(0,'dir','horz','linewidth',lw,'linestyle','--'); drawline(12.5,'dir','vert','linewidth',28,'color',[1 1 1]); drawline(12.5,'dir','vert','linewidth',lw);
            xlabel('Block number'); ylabel('Mean Delta MT (ms)'); set(gca,'fontsize',fs);
            % error learning curves
            D5=tapply(D,{'SN','BN','repNum'},...
                {D.isError,'sum','name','sumErr'},...
                {D.TN,'numel','name','nTrials'},...
                'subset',D.dummy==0);
            D5.ER=((D5.sumErr./D5.nTrials)*100);
            subplot(3,4,[7,8]); title('Learning: Error rate (ER)');
            lineplot(D5.BN,D5.ER,'split',D5.repNum,'errorbars','plusminus','errorcap',0,'linewidth',lw,'linecolor',cr,'leg',lr,'leglocation','northeast'); drawline(15,'dir','horz','linewidth',lw,'linestyle','--');
            xlabel('Block number'); ylabel('Mean ER (%)'); set(gca,'fontsize',fs); set(gca,'ylim',[0 100]); drawline(12.5,'dir','vert','linewidth',28,'color',[1 1 1]); drawline(12.5,'dir','vert','linewidth',lw);
            D6=tapply(D5,{'SN','BN'},...
                {D5.ER,'mean','name','ERfirst','subset',D5.repNum==1},...
                {D5.ER,'mean','name','ERrest','subset',D5.repNum>1});
            subplot(3,4,[11,12]); title('Delta ER Rep1-RepOther');
            lineplot(D6.BN,D6.ERfirst-D6.ERrest,'errorbars','shade','shadecolor',[.5 .5 .5],'linewidth',lw,'linecolor',[.5 .5 .5],'leg',lr,'leglocation','northeast'); drawline(0,'dir','horz','linewidth',lw,'linestyle','--');
            xlabel('Block number'); ylabel('Mean Delta ER (%)'); set(gca,'fontsize',fs); drawline(12.5,'dir','vert','linewidth',28,'color',[1 1 1]); drawline(12.5,'dir','vert','linewidth',lw);
            % scatter plots MT: relationship bwetween MT and repetition effect
            figure('Name',sprintf('Scatter plots MT - group (N=%d)',numel(subj)));
            D7=tapply(D,{'SN'},...
                {D.MT,'mean','name','mMT'},...
                {D.MT,'mean','name','MTfirst','subset',D.repNum==1},...
                {D.MT,'mean','name','MTrest','subset',D.repNum>1},...
                'subset',D.dummy==0 & D.isError==0);
            subplot(2,2,1); title('All subjects, both days'); hold on;
            scatterplot(D7.mMT,D7.MTfirst-D7.MTrest,'regression','linear',...
                'markersize',9,'label',subj,'printcorr');
            xlabel('Mean movement time (MT)'); ylabel('Repetition effect (Delta MT)'); set(gca,'fontsize',fs); axis([min(D7.mMT)-100 max(D7.mMT)+100 min(D7.MTfirst-D7.MTrest)-20 max(D7.MTfirst-D7.MTrest)+20]);
            D7=tapply(D,{'SN','day'},...
                {D.MT,'mean','name','mMT'},...
                {D.MT,'mean','name','MTfirst','subset',D.repNum==1},...
                {D.MT,'mean','name','MTrest','subset',D.repNum>1},...
                'subset',D.dummy==0 & D.isError==0);
            subplot(2,2,3); title('All subjects, each day'); hold on;
            scatterplot(D7.mMT,D7.MTfirst-D7.MTrest,'split',D7.day,'regression','linear',...
                'markertype','o','markercolor',cd,'markerfill',cd,'markersize',9,'leg',dayleg,'leglocation','northeast',...
                'label',subl,'printcorr');
            xlabel('Mean movement time (MT)'); ylabel('Repetition effect (Delta MT)'); set(gca,'fontsize',fs); axis([min(D7.mMT)-100 max(D7.mMT)+100 min(D7.MTfirst-D7.MTrest)-20 max(D7.MTfirst-D7.MTrest)+20]);
            D8=tapply(D,{'seqNum'},...
                {D.MT,'mean','name','mMT'},...
                {D.MT,'mean','name','MTfirst','subset',D.repNum==1},...
                {D.MT,'mean','name','MTrest','subset',D.repNum>1},...
                'subset',D.dummy==0 & D.isError==0);
            subplot(2,2,2); title('All sequences, both days'); hold on;
            scatterplot(D8.mMT,D8.MTfirst-D8.MTrest,'regression','linear',...
                'markersize',9,'label',1:8,'printcorr');
            xlabel('Mean movement time (MT)'); ylabel('Repetition effect (Delta MT)'); set(gca,'fontsize',fs); axis([min(D8.mMT)-50 max(D8.mMT)+50 min(D8.MTfirst-D8.MTrest)-10 max(D8.MTfirst-D8.MTrest)+10]);
            D8=tapply(D,{'seqNum','day'},...
                {D.MT,'mean','name','mMT'},...
                {D.MT,'mean','name','MTfirst','subset',D.repNum==1},...
                {D.MT,'mean','name','MTrest','subset',D.repNum>1},...
                'subset',D.dummy==0 & D.isError==0);
            subplot(2,2,4); title('All sequences, each day'); hold on;
            scatterplot(D8.mMT,D8.MTfirst-D8.MTrest,'split',D8.day,'regression','linear',...
                'markertype','o','markercolor',cd,'markerfill',cd,'markersize',9,'leg',dayleg,'leglocation','northeast',...
                'label',seql,'printcorr');
            xlabel('Mean movement time (MT)'); ylabel('Repetition effect (Delta MT)'); set(gca,'fontsize',fs); axis([min(D8.mMT)-50 max(D8.mMT)+50 min(D8.MTfirst-D8.MTrest)-10 max(D8.MTfirst-D8.MTrest)+10]);
            % scatter plots RT: relationship bwetween RT and repetition effect
            figure('Name',sprintf('Scatter plots RT - group (N=%d)',numel(subj)));
            D9=tapply(D,{'SN','day'},...
                {D.RT,'mean','name','mRT'},...
                {D.RT,'mean','name','RTfirst','subset',D.repNum==1},...
                {D.RT,'mean','name','RTrest','subset',D.repNum>1},...
                'subset',D.dummy==0 & D.isError==0);
            subplot(2,2,1); title('All subjects, each day'); hold on;
            scatterplot(D9.mRT,D9.RTfirst-D9.RTrest,'split',D9.day,'regression','linear',...
                'markertype','o','markercolor',cd,'markerfill',cd,'markersize',9,'leg',dayleg,'leglocation','northeast',...
                'label',subl,'printcorr');
            xlabel('Mean reaction time (RT)'); ylabel('Repetition effect (Delta RT)'); set(gca,'fontsize',fs); axis([min(D9.mRT)-50 max(D9.mRT)+50 min(D9.RTfirst-D9.RTrest)-20 max(D9.RTfirst-D9.RTrest)+20]);
            D9=tapply(D,{'SN','day'},...
                {D.MT,'mean','name','mMT'},...
                {D.RT,'mean','name','mRT'},...
                'subset',D.dummy==0 & D.isError==0);
            subplot(2,2,2); title('All subjects, each day'); hold on;
            scatterplot(D9.mRT,D9.mMT,'split',D9.day,'regression','linear',...
                'markertype','o','markercolor',cd,'markerfill',cd,'markersize',9,'leg',dayleg,'leglocation','northeast',...
                'label',subl,'printcorr');
            xlabel('Mean reaction time (RT)'); ylabel('Mean movement time (MT)'); set(gca,'fontsize',fs); axis([min(D9.mRT)-50 max(D9.mRT)+50 min(D9.mMT)-100 max(D9.mMT)+100]);
            D9=tapply(D,{'SN','BN'},...
                {D.MT,'mean','name','mMT'},...
                {D.RT,'mean','name','mRT'},...
                'subset',D.dummy==0 & D.isError==0);
            subplot(2,2,[3,4]); title('Learning: Movement time (MT) vs Reaction time (RT)');
            lineplot(D9.BN,D9.mMT,'errorbars','shade','shadecolor',[.3 .3 .3],'linewidth',lw,'linecolor',[.3 .3 .3],'leg',{'MT'}); hold on;
            lineplot(D9.BN,D9.mRT,'errorbars','shade','shadecolor',[.7 .7 .7],'linewidth',lw,'linecolor',[.7 .7 .7],'leg',{'RT'}); hold on;
            drawline(12.5,'dir','vert','linewidth',28,'color',[1 1 1]); drawline(12.5,'dir','vert','linewidth',lw);
            xlabel('Block number'); ylabel('Mean MT | RT (ms)'); set(gca,'fontsize',fs);
            
            % out
            D.D1=D1; D.D2=D2; D.D3=D3; D.D4=D4; D.D5=D5; D.D6=D6; D.D7=D7; D.D8=D8; D.D9=D9; %incorporate the sub-structures as fields of main structure
            varargout={D}; %return main structure
        end
        
    case 'pointsAnalysis' %proportion of points/errors (effect of thresholds) per block number,subject and group
        c={[1 0 0],[.9 .9 .9],[.6 .6 .6],[.1 .1 .1]};
        l={'Err','0','+1','+3'};
        allD=[];
        figure;
        for s=1:length(subj)
            D=load(fullfile(pathToData,sprintf('sr1_%s.mat',subj{s})));
            D.points(D.isError==1)=-1;
            D=tapply(D,{'BN','points'},...
                {D.points,'numel','name','nPoints'},...
                'subset',D.dummy==0);
            subplot(ceil(length(subj)/2),2,s);
            lineplot(D.BN,D.nPoints,'split',D.points,'linewidth',lw,'linecolor',c,'leg',l,'leglocation','northeast');
            title(sprintf('Subject %s',subj{s})); xlabel('BN'); ylabel('N trials'); set(gca,'fontsize',fs); drawline(12.5,'dir','vert','linewidth',28,'color',[1 1 1]); drawline(12.5,'dir','vert','linewidth',lw); set(gca,'ylim',[0 50]);
            allD=addstruct(allD,D);
        end
        figure;
        lineplot(allD.BN,allD.nPoints,'split',allD.points,'errorbars','shade','shadecolor',c,'linewidth',lw,'linecolor',c,'leg',l,'leglocation','northeast');
        title('Group average'); xlabel('BN'); ylabel('N trials per block'); set(gca,'fontsize',fs); drawline(12.5,'dir','vert','linewidth',64,'color',[1 1 1]); drawline(12.5,'dir','vert','linewidth',lw);
        figure;
        for s=1:length(subj)
            D=load(fullfile(pathToData,sprintf('sr1_%s.mat',subj{s})));
            D.points(D.isError==1)=-1;
            D=tapply(D,{'points'},...
                {D.points,'numel','name','nPoints'},...
                'subset',D.dummy==0);
            subplot(2,ceil(length(subj)/2),s);
            barplot((1:numel(unique(D.points)))',D.nPoints/sum(D.nPoints),'split',D.points,'facecolor',c); xticklabels({'Err','0','+1','+3'});
            title(sprintf('Subject %s',subj{s})); xlabel('Points'); ylabel('Proportion'); set(gca,'fontsize',fs); set(gca,'ylim',[0 .8]);
        end
        
    case 'repEff_ET' % analysis of repetition effect on ET
        if nargin>1 % load single subj data
            subj = varargin{1};
            D = load( fullfile(pathToData, sprintf('sr1_%s.mat', subj)));
        else % load group data
            D = load( fullfile(pathToAnalyze, 'sr1_all_data.mat'));
        end
        
        %-------------------------------------------------------------------------------------------------------------------------------------
        % select repetitions of interest
        if numel(rs)>1; D = getrow(D, ismember(D.repNum, rs)); end
        
        %-------------------------------------------------------------------------------------------------------------------------------------
        % open figure
        if nargin>1; figure('Name',sprintf('Repetition effect - subj %02d',str2double(varargin{1}(2:3)))); else; figure('Name',sprintf('Repetition effect - group (N=%d)',ns)); end
        set(gcf, 'Units','normalized', 'Position',[0.1,0.1,0.8,0.8], 'Resize','off', 'Renderer','painters');
        
        %-------------------------------------------------------------------------------------------------------------------------------------
        % repEff create summary table
        T = tapply(D, {'SN', 'day', 'isRep'}, ...
            {D.ET, 'nanmedian', 'name', 'ET'}, ...
            'subset', D.isError==0);
        
        %-------------------------------------------------------------------------------------------------------------------------------------
        % normalize MT data to remove between-subject variability (i.e. plot within-subject standard error)
        T = normData(T, {'ET'}, 'sub');
        
        %-------------------------------------------------------------------------------------------------------------------------------------
        % make sure that you have one value per subject for each condition
        %pivottable(T.isRep, T.SN, T.normET, 'length');
        
        %-------------------------------------------------------------------------------------------------------------------------------------
        % repEff box plot
        subplot(2,2,1); title('Repetition Effect'); hold on;
        plt.box(T.day, T.normET, 'split',T.isRep, 'style',isrepstybox, 'leg',isrepleg, 'leglocation','northeast');
        
        %-------------------------------------------------------------------------------------------------------------------------------------
        % repEff line plot
        T.isRep(T.isRep == 1) = 2; T.isRep(T.isRep == 0) = 1;
        hold on;
        plt.line([T.day T.isRep], T.normET, 'split', T.SN, 'errorbars','none', 'style',lightgraysty, 'leg','skip');
        hold on;
        plt.line([T.day T.isRep], T.normET, 'errorbars','plusminus', 'style',blacksty, 'leg','skip');
        xt = xticks; xticks([sum(xt(1:2))/2 sum(xt(3:4))/2]); xticklabels(dayleg);
        ylabel('Execution time (ms)'); set(gca,'fontsize',fs); axis square;
        
        %-------------------------------------------------------------------------------------------------------------------------------------
        % stats
        ttest(T.ET(T.isRep==1 & T.day==1), T.ET(T.isRep==2 & T.day==1), 2, 'paired');
        ttest(T.ET(T.isRep==1 & T.day==2), T.ET(T.isRep==2 & T.day==2), 2, 'paired');
        
        % -------------------------------------------------------------------------------------------------------------------------------------
        % Rep num create summary table
        if maxRep > 0; D.repNum(D.repNum >= maxRep) = maxRep; end % put a ceiling to nReps
        T = tapply(D, {'SN', 'repNum'}, ...
            {D.ET, 'nanmedian', 'name', 'ET'}, ...
            'subset', D.isError==0);
        
        %-------------------------------------------------------------------------------------------------------------------------------------
        % normalize MT data to remove between-subject variability (i.e. plot within-subject standard error)
        T = normData(T, {'ET'}, 'sub');
        
        %-------------------------------------------------------------------------------------------------------------------------------------
        % make sure that you have one value per subject for each condition
        %pivottable(T.repNum, T.SN, T.normET, 'length');
        
        %-------------------------------------------------------------------------------------------------------------------------------------
        % Rep num line plot
        subplot(2,2,2); title('Repetition number');
        [~,~] = plt.line(T.repNum, T.normET, 'errorbars','shade', 'style',darkgraysty, 'leg','skip',  'leglocation','northeast');
        xlabel('Repetition number'); ylabel('ET (ms)'); set(gca,'fontsize',fs); axis square;
        %ylim([690 810]);
        labels = cell(1, max(D.repNum)+1); labels{1, 1} = 'Switch';
        for l = 1 : max(D.repNum); labels{1, l+1} = num2str(l); end
        labels{1, end} = sprintf('%s', labels{1, end}); xticklabels(labels);
        
        %-------------------------------------------------------------------------------------------------------------------------------------
        nq = 11; % number of quantiles
        [D.I,D.M,~] = split_data(D.ET, 'split',[D.SN D.isRep], 'numquant',nq, 'subset',D.isError==0);
        T = tapply(D, {'SN', 'I', 'isRep'}, ...
            {D.ET, 'nanmedian', 'name', 'ET'}, ...
            'subset',D.isError==0);
        
        % normalize MT data to remove between-subject variability (i.e. plot within-subject standard error)
        T = normData(T, {'ET'}, 'sub');
        
        subplot(2,2,3); title(''); hold on;
        plt.line(T.I, T.normET, 'split',T.isRep, 'style',isrepsty, 'leg',isrepleg, 'leglocation','northwest');
        xticklabels(linspace(0,100,nq));
        xlabel('ET percentile (%)'); ylabel('ET (ms)');  set(gca,'fontsize',fs); axis square;
        xlim([0.5 nq+0.5]);
        ylim([380 1140]);
        drawline(round(nq/2), 'dir','vert', 'linestyle','--');
        
        %-------------------------------------------------------------------------------------------------------------------------------------
        [D.I,D.M,~] = split_data(D.ET, 'split',[D.SN D.isRep], 'numquant',nq, 'subset',D.isError==0);
        T = tapply(D, {'SN', 'I'}, ...
            {D.ET, 'nanmedian', 'name', 'ET'}, ...
            {D.ET, 'nanmedian', 'name', 'ETs', 'subset',D.isRep==0}, ...
            {D.ET, 'nanmedian', 'name', 'ETr', 'subset',D.isRep==1}, ...
            'subset',D.isError==0);
        
        % normalize MT data to remove between-subject variability (i.e. plot within-subject standard error)
        T = normData(T, {'ETs', 'ETr', 'ET'}, 'sub');
        
        subplot(2,2,4); title(''); hold on;
        plt.line(T.I, (nanplus(T.normETs,-T.normETr)./T.normET)*100, 'style',darkgraysty);
        xticklabels(linspace(0,100,nq));
        xlabel('ET percentile (%)'); ylabel('Repetition difference (% of ET)');  set(gca,'fontsize',fs); axis square;
        xlim([0.5 nq+0.5]);
        drawline(round(nq/2), 'dir','vert', 'linestyle','--');
        
        %-------------------------------------------------------------------------------------------------------------------------------------
        % out
        D.T=T; %incorporate the sub-structures as fields of main structure
        varargout={D}; %return main structure
        
    case 'repEff_RT' % analysis of repetition effect on RT
        if nargin>1 % load single subj data
            subj = varargin{1};
            D = load( fullfile(pathToData, sprintf('sr1_%s.mat', subj)));
        else % load group data
            D = load( fullfile(pathToAnalyze, 'sr1_all_data.mat'));
        end
        
        %-------------------------------------------------------------------------------------------------------------------------------------
        % select repetitions of interest
        if numel(rs)>1; D = getrow(D, ismember(D.repNum, rs)); end
        
        %-------------------------------------------------------------------------------------------------------------------------------------
        % open figure
        if nargin>1; figure('Name',sprintf('Repetition effect - subj %02d',str2double(varargin{1}(2:3)))); else; figure('Name',sprintf('Repetition effect - group (N=%d)',ns)); end
        set(gcf, 'Units','normalized', 'Position',[0.1,0.1,0.8,0.8], 'Resize','off', 'Renderer','painters');
        
        %-------------------------------------------------------------------------------------------------------------------------------------
        % repEff create summary table
        T = tapply(D, {'SN', 'day', 'isRep'}, ...
            {D.RT, 'nanmedian', 'name', 'RT'}, ...
            'subset', D.isError==0);
        
        %-------------------------------------------------------------------------------------------------------------------------------------
        % normalize MT data to remove between-subject variability (i.e. plot within-subject standard error)
        T = normData(T, {'RT'}, 'sub');
        
        %-------------------------------------------------------------------------------------------------------------------------------------
        % make sure that you have one value per subject for each condition
        %pivottable(T.isRep, T.SN, T.normRT, 'length');
        
        %-------------------------------------------------------------------------------------------------------------------------------------
        % repEff box plot
        subplot(2,2,1); title('Repetition Effect'); hold on;
        plt.box(T.day, T.normRT, 'split',T.isRep, 'style',isrepstybox, 'leg',isrepleg, 'leglocation','northeast');
        
        %-------------------------------------------------------------------------------------------------------------------------------------
        % repEff line plot
        T.isRep(T.isRep == 1) = 2; T.isRep(T.isRep == 0) = 1;
        hold on;
        plt.line([T.day T.isRep], T.normRT, 'split', T.SN, 'errorbars','none', 'style',lightgraysty, 'leg','skip');
        hold on;
        plt.line([T.day T.isRep], T.normRT, 'errorbars','plusminus', 'style',blacksty, 'leg','skip');
        xt = xticks; xticks([sum(xt(1:2))/2 sum(xt(3:4))/2]); xticklabels(dayleg);
        ylabel('Reaction time (ms)'); set(gca,'fontsize',fs); axis square;
        
        %-------------------------------------------------------------------------------------------------------------------------------------
        % stats
        ttest(T.RT(T.isRep==1 & T.day==1), T.RT(T.isRep==2 & T.day==1), 2, 'paired');
        ttest(T.RT(T.isRep==1 & T.day==2), T.RT(T.isRep==2 & T.day==2), 2, 'paired');
        
        % -------------------------------------------------------------------------------------------------------------------------------------
        % Rep num create summary table
        if maxRep > 0; D.repNum(D.repNum >= maxRep) = maxRep; end % put a ceiling to nReps
        T = tapply(D, {'SN', 'repNum'}, ...
            {D.RT, 'nanmedian', 'name', 'RT'}, ...
            'subset', D.isError==0);
        
        %-------------------------------------------------------------------------------------------------------------------------------------
        % normalize MT data to remove between-subject variability (i.e. plot within-subject standard error)
        T = normData(T, {'RT'}, 'sub');
        
        %-------------------------------------------------------------------------------------------------------------------------------------
        % make sure that you have one value per subject for each condition
        %pivottable(T.repNum, T.SN, T.normRT, 'length');
        
        %-------------------------------------------------------------------------------------------------------------------------------------
        % Rep num line plot
        subplot(2,2,2); title('Repetition number');
        [~,~] = plt.line(T.repNum, T.normRT, 'errorbars','shade', 'style',darkgraysty, 'leg','skip',  'leglocation','northeast');
        xlabel('Repetition number'); ylabel('RT (ms)'); set(gca,'fontsize',fs); axis square;
        %ylim([690 810]);
        labels = cell(1, max(D.repNum)+1); labels{1, 1} = 'Switch';
        for l = 1 : max(D.repNum); labels{1, l+1} = num2str(l); end
        labels{1, end} = sprintf('%s', labels{1, end}); xticklabels(labels);
        
        %-------------------------------------------------------------------------------------------------------------------------------------
        nq = 11; % number of quantiles
        [D.I,D.M,~] = split_data(D.ET, 'split',[D.SN D.isRep], 'numquant',nq, 'subset',D.isError==0);
        T = tapply(D, {'SN', 'I', 'isRep'}, ...
            {D.RT, 'nanmedian', 'name', 'RT'}, ...
            'subset',D.isError==0);
        
        % normalize MT data to remove between-subject variability (i.e. plot within-subject standard error)
        T = normData(T, {'RT'}, 'sub');
        
        subplot(2,2,3); title(''); hold on;
        plt.line(T.I, T.normRT, 'split',T.isRep, 'style',isrepsty, 'leg',isrepleg, 'leglocation','northwest');
        xticklabels(linspace(0,100,nq));
        xlabel('RT percentile (%)'); ylabel('RT (ms)');  set(gca,'fontsize',fs); axis square;
        xlim([0.5 nq+0.5]);
        ylim([351 551]);
        drawline(round(nq/2), 'dir','vert', 'linestyle','--');
        
        %-------------------------------------------------------------------------------------------------------------------------------------
        [D.I,D.M,~] = split_data(D.RT, 'split',[D.SN D.isRep], 'numquant',nq, 'subset',D.isError==0);
        T = tapply(D, {'SN', 'I'}, ...
            {D.RT, 'nanmedian', 'name', 'RT'}, ...
            {D.RT, 'nanmedian', 'name', 'RTs', 'subset',D.isRep==0}, ...
            {D.RT, 'nanmedian', 'name', 'RTr', 'subset',D.isRep==1}, ...
            'subset',D.isError==0);
        
        % normalize MT data to remove between-subject variability (i.e. plot within-subject standard error)
        T = normData(T, {'RTs', 'RTr', 'RT'}, 'sub');
        
        subplot(2,2,4); title(''); hold on;
        plt.line(T.I, (nanplus(T.normRTs,-T.normRTr)./T.normRT)*100, 'style',darkgraysty);
        xticklabels(linspace(0,100,nq));
        xlabel('RT percentile (%)'); ylabel('Repetition difference (% of RT)');  set(gca,'fontsize',fs); axis square;
        xlim([0.5 nq+0.5]);
        drawline(round(nq/2), 'dir','vert', 'linestyle','--');
        
        %-------------------------------------------------------------------------------------------------------------------------------------
        % out
        D.T=T; %incorporate the sub-structures as fields of main structure
        varargout={D}; %return main structure
        
    case 'repDiff' % analysis of Swc-Rep difference on both ET and RT
        if nargin>1 % load single subj data
            subj = varargin{1};
            D = load( fullfile(pathToData, sprintf('sr1_%s.mat', subj)));
        else % load group data
            D = load( fullfile(pathToAnalyze, 'sr1_all_data.mat'));
        end
        
        %-------------------------------------------------------------------------------------------------------------------------------------
        % select repetitions of interest
        if numel(rs)>1; D = getrow(D, ismember(D.repNum, rs)); end
        
        %-------------------------------------------------------------------------------------------------------------------------------------
%         % open figure
%         if nargin>1; figure('Name',sprintf('Repetition effect - subj %02d',str2double(varargin{1}(2:3)))); else; figure('Name',sprintf('Repetition effect - group (N=%d)',ns)); end
%         set(gcf, 'Units','normalized', 'Position',[0.1,0.1,0.8,0.8], 'Resize','off', 'Renderer','painters');
        
        %-------------------------------------------------------------------------------------------------------------------------------------
        % repDiff create summary table
        T = tapply(D, {'SN', 'BN', 'day', 'isRep'}, ...
            {D.ET, 'nanmedian', 'name', 'ET'}, ...
            {D.RT, 'nanmedian', 'name', 'RT'}, ...
            'subset', D.isError==0);
        
        %-------------------------------------------------------------------------------------------------------------------------------------
        % normalize MT data to remove between-subject variability (i.e. plot within-subject standard error)
        T = normData(T, {'ET', 'RT'}, 'sub');
        
        %-------------------------------------------------------------------------------------------------------------------------------------
        % repEff line plot
        subplot(2,2,3); title(''); hold on;
        plt.line([T.day T.BN], T.ET, 'split',[T.isRep], 'style',isrepsty, 'leg',isrepleg, 'leglocation','northeast');
        xlabel('Block number'); ylabel('Execution time (ms)'); set(gca,'fontsize',fs); axis square;
        %xt = xticks; xticks([xt(6), xt(18)]); xticklabels(dayleg); 
        ylim([480 1020]);
        %drawline(9.2, 'dir','vert', 'linestyle',':');
        drawline(12.5, 'dir','vert', 'linestyle',':');
        
        % stats
        ttest(T.ET(T.BN==1 & T.isRep==0), T.ET(T.BN==1 & T.isRep==1), 2, 'paired');
        ttest(T.ET(T.BN==2 & T.isRep==0), T.ET(T.BN==2 & T.isRep==1), 2, 'paired');
        ttest(T.ET(T.BN==3 & T.isRep==0), T.ET(T.BN==3 & T.isRep==1), 2, 'paired');
        ttest(T.ET(T.BN==4 & T.isRep==0), T.ET(T.BN==4 & T.isRep==1), 2, 'paired');
        
        ttest(T.ET(T.BN==13 & T.isRep==0), T.ET(T.BN==13 & T.isRep==1), 2, 'paired');
        ttest(T.ET(T.BN==14 & T.isRep==0), T.ET(T.BN==14 & T.isRep==1), 2, 'paired');
        ttest(T.ET(T.BN==15 & T.isRep==0), T.ET(T.BN==15 & T.isRep==1), 2, 'paired');
        ttest(T.ET(T.BN==16 & T.isRep==0), T.ET(T.BN==16 & T.isRep==1), 2, 'paired');
        
%         subplot(2,2,3); title(''); hold on;
%         plt.line([T.day T.BN], T.normRT, 'split',[T.isRep], 'style',isrepsty, 'leg',isrepleg, 'leglocation','northeast');
%         ylabel('Reaction time (ms)'); set(gca,'fontsize',fs); axis square;
%         xt = xticks; xticks([xt(6), xt(18)]); xticklabels(dayleg);  ylim([340 610]);
%         drawline(9.2, 'dir','vert', 'linestyle',':');
        
        %-------------------------------------------------------------------------------------------------------------------------------------
        % repDiff create summary table
        T = tapply(D, {'SN', 'BN', 'day'}, ...
            {D.ET, 'nanmedian', 'name', 'ET'}, ...
            {D.ET, 'nanmedian', 'name', 'ETswc', 'subset',D.isRep==0}, ...
            {D.ET, 'nanmedian', 'name', 'ETrep', 'subset',D.isRep==1}, ...
            {D.RT, 'nanmedian', 'name', 'RT'}, ...
            {D.RT, 'nanmedian', 'name', 'RTswc', 'subset',D.isRep==0}, ...
            {D.RT, 'nanmedian', 'name', 'RTrep', 'subset',D.isRep==1}, ...
            'subset', D.isError==0);
        
        %-------------------------------------------------------------------------------------------------------------------------------------
        % normalize MT data to remove between-subject variability (i.e. plot within-subject standard error)
        T = normData(T, {'ET', 'ETswc', 'ETrep', 'RT', 'RTswc', 'RTrep'}, 'sub');
        
        %-------------------------------------------------------------------------------------------------------------------------------------
        % Diff line plot
        subplot(2,2,4);
        %plt.scatter([T.BN], ((T.ETswc-T.ETrep)./T.ET)*100, 'split',T.day, 'leg','skip', 'leglocation','northeast');
        %plt.scatter([T.BN], T.ETswc-T.ETrep, 'split',T.day, 'style',isrepsty, 'leg','skip', 'leglocation','northeast');
        plt.line([T.day T.BN], ((T.ETswc-T.ETrep)./T.ET)*100, 'style',graysty, 'leg','skip', 'leglocation','northeast');
        hold on;
        plt.scatter([T.BN], ((T.ETswc-T.ETrep)./T.ET)*100, 'split',T.day, 'style',blacksty, 'leg','skip', 'leglocation','northeast');
        %plt.line([T.day T.BN], (T.ETswc-T.ETrep), 'style',graysty, 'leg','skip', 'leglocation','northeast');
        xlabel('Block number'); ylabel('Repetition difference (% of ET)'); set(gca,'fontsize',fs); axis square;
        %xt = xticks; xticks([xt(6), xt(18)]); xticklabels(dayleg); 
        ylim([-2 8]); %ylim([-15 65]);
        drawline(0, 'dir','horz', 'linestyle','--'); %drawline(9.2, 'dir','vert', 'linestyle',':');
        drawline(12.5, 'dir','vert', 'linestyle',':');
        
        T = tapply(D, {'SN', 'BN'}, ...
            {D.ET, 'nanmedian', 'name', 'ET'}, ...
            {D.ET, 'nanmedian', 'name', 'ETswc', 'subset',D.isRep==0}, ...
            {D.ET, 'nanmedian', 'name', 'ETrep', 'subset',D.isRep==1}, ...
            'subset', D.isError==0);
        
        % stats
        ET = ((T.ETswc-T.ETrep)./T.ET)*100; 
        ttest(ET(T.BN==1), 0, 2, 'onesample');
        ttest(ET(T.BN==2), 0, 2, 'onesample');
        ttest(ET(T.BN==3), 0, 2, 'onesample');
        ttest(ET(T.BN==4), 0, 2, 'onesample');
        
        ttest(ET(T.BN==13), 0, 2, 'onesample');
        ttest(ET(T.BN==14), 0, 2, 'onesample');
        ttest(ET(T.BN==15), 0, 2, 'onesample');
        ttest(ET(T.BN==16), 0, 2, 'onesample');
        
%         subplot(2,2,4);
%         plt.line([T.day T.BN], ((T.normRTswc-T.normRTrep)./T.normET)*100, 'split',[], 'style',graysty, 'leg','skip', 'leglocation','northeast');
%         ylabel('Repetition benefit (% of RT)'); set(gca,'fontsize',fs); axis square;
%         xt = xticks; xticks([xt(6), xt(18)]); xticklabels(dayleg); %ylim([-20 160]);
%         drawline(0, 'dir','horz', 'linestyle','--', 'linewidth',lw); drawline(9.2, 'dir','vert', 'linestyle',':');
        
        %-------------------------------------------------------------------------------------------------------------------------------------
        % stats
        %ttest(T.RT(T.isRep==0), T.RT(T.isRep==1), 2, 'paired');
        
        %-------------------------------------------------------------------------------------------------------------------------------------
        % out
        D.T=T; %incorporate the sub-structures as fields of main structure
        varargout={D}; %return main structure
        
    case 'repCorr' % analysis of repEff correlation with either ET or RT
        if nargin>1 % load single subj data
            subj = varargin{1};
            D = load( fullfile(pathToData, sprintf('sr1_%s.mat', subj)));
        else % load group data
            D = load( fullfile(pathToAnalyze, 'sr1_all_data.mat'));
        end
        
        %-------------------------------------------------------------------------------------------------------------------------------------
        % select repetitions of interest
        if numel(rs)>1; D = getrow(D, ismember(D.repNum, rs)); end
        
        %-------------------------------------------------------------------------------------------------------------------------------------
        % open figure
        if nargin>1; figure('Name',sprintf('Repetition effect - subj %02d',str2double(varargin{1}(2:3)))); else; figure('Name',sprintf('Repetition effect - group (N=%d)',ns)); end
        set(gcf, 'Units','normalized', 'Position',[0.1,0.1,0.8,0.8], 'Resize','off', 'Renderer','painters');
        
        %-------------------------------------------------------------------------------------------------------------------------------------
        % repDiff create summary table
        T = tapply(D, {'SN', 'day'}, ...
            {D.ET, 'nanmedian', 'name', 'ETswc', 'subset',D.isRep==0}, ...
            {D.ET, 'nanmedian', 'name', 'ETrep', 'subset',D.isRep==1}, ...
            {D.RT, 'nanmedian', 'name', 'RTswc', 'subset',D.isRep==0}, ...
            {D.RT, 'nanmedian', 'name', 'RTrep', 'subset',D.isRep==1}, ...
            {D.ET, 'nanmedian', 'name', 'ET'}, ...
            {D.RT, 'nanmedian', 'name', 'RT'}, ...
            'subset', D.isError==0);
        
        %-------------------------------------------------------------------------------------------------------------------------------------
        % normalize MT data to remove between-subject variability (i.e. plot within-subject standard error)
        T = normData(T, {'ETswc', 'ETrep', 'RTswc', 'RTrep', 'ET', 'RT'}, 'sub');
        
        %-------------------------------------------------------------------------------------------------------------------------------------
        % repEff scatter plot
        subplot(2,2,3); title('Execution time'); hold on;
        plt.scatter(T.ET, T.normETswc-T.normETrep, 'split',[T.day], 'style',daysty, 'label',reshape(repmat(subj,2,1),1,ns*2), 'printcorr',1, 'leg',dayleg);%, 'subset',T.day == 1);
        xlabel('ET'); ylabel('Repetition effect on ET (ms)'); set(gca,'fontsize',fs); axis square;
        xlim([100 1700]); ylim([-30 55]);
        
        subplot(2,2,4); title('Reaction time'); hold on;
        plt.scatter(T.RT, T.normRTswc-T.normRTrep, 'split',[T.day], 'style',daysty, 'label',reshape(repmat(subj,2,1),1,ns*2), 'printcorr',1, 'leg',dayleg);%, 'subset',T.day == 1);
        xlabel('RT'); ylabel('Repetition effect on RT (ms)'); set(gca,'fontsize',fs); axis square;
        xlim([200 800]); ylim([-100 200]);
        
        %-------------------------------------------------------------------------------------------------------------------------------------
        % stats
        %ttest(T.RT(T.isRep==0), T.RT(T.isRep==1), 2, 'paired');
        
        %-------------------------------------------------------------------------------------------------------------------------------------
        % out
        D.T=T; %incorporate the sub-structures as fields of main structure
        varargout={D}; %return main structure
        
    case 'repSff' % analysis of repEff interaction with same first finger on ET and RT
        if nargin>1 % load single subj data
            subj = varargin{1};
            D = load( fullfile(pathToData, sprintf('sr1_%s.mat', subj)));
        else % load group data
            D = load( fullfile(pathToAnalyze, 'sr1_all_data.mat'));
        end
        
        %-------------------------------------------------------------------------------------------------------------------------------------
        % select repetitions of interest
        if numel(rs)>1; D = getrow(D, ismember(D.repNum, rs)); end
        
        %-------------------------------------------------------------------------------------------------------------------------------------
        % open figure
        if nargin>1; figure('Name',sprintf('Repetition effect - subj %02d',str2double(varargin{1}(2:3)))); else; figure('Name',sprintf('Repetition effect - group (N=%d)',ns)); end
        set(gcf, 'Units','normalized', 'Position',[0.1,0.1,0.8,0.8], 'Resize','off', 'Renderer','painters');
        
        %-------------------------------------------------------------------------------------------------------------------------------------
        % repDiff create summary table
        T = tapply(D, {'SN', 'day', 'sff'}, ...
            {D.ET, 'nanmedian', 'name', 'ET'}, ...
            {D.RT, 'nanmedian', 'name', 'RT'}, ...
            'subset', D.isError==0);
        
        %-------------------------------------------------------------------------------------------------------------------------------------
        % normalize MT data to remove between-subject variability (i.e. plot within-subject standard error)
        T = normData(T, {'ET', 'RT'}, 'sub');
        
        %-------------------------------------------------------------------------------------------------------------------------------------
        % repEff scatter plot
        subplot(2,2,1); title('Day 1'); hold on;
        plt.bar(T.sff, T.normET, 'split',[], 'style',darkgraysty, 'subset',T.day == 1);
        xlabel(''); xticklabels({'Swc', 'Sff', 'Rep'}); ylabel('Execution time (ms)'); set(gca,'fontsize',fs); axis square;
        ylim([500 800]);
        subplot(2,2,2); title('Day 2'); hold on;
        plt.bar(T.sff, T.normET, 'split',[], 'style',darkgraysty, 'subset',T.day == 2);
        xlabel(''); xticklabels({'Swc', 'Sff', 'Rep'}); ylabel('ET (ms)'); set(gca,'fontsize',fs); axis square;
        ylim([500 800]);
        
        subplot(2,2,3); title('Day 1'); hold on;
        plt.bar(T.sff, T.normRT, 'split',[], 'style',graysty, 'subset',T.day == 1);
        xlabel(''); xticklabels({'Swc', 'Sff', 'Rep'}); ylabel('Reaction time (ms)'); set(gca,'fontsize',fs); axis square;
        ylim([300 600]);
        subplot(2,2,4); title('Day 2'); hold on;
        plt.bar(T.sff, T.normRT, 'split',[], 'style',graysty, 'subset',T.day == 2);
        xlabel(''); xticklabels({'Swc', 'Sff', 'Rep'}); ylabel('RT (ms)'); set(gca,'fontsize',fs); axis square;
        ylim([300 600]);
        
        %-------------------------------------------------------------------------------------------------------------------------------------
        % stats
        %ttest(T.RT(T.isRep==0), T.RT(T.isRep==1), 2, 'paired');
        
        %-------------------------------------------------------------------------------------------------------------------------------------
        % out
        D.T=T; %incorporate the sub-structures as fields of main structure
        varargout={D}; %return main structure
        
    otherwise
        error('no such case!')
end;