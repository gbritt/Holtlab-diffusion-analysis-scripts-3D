function [traceList,info]=tracksFinal_2_traceList(tracksFinal, track_pars)
%
% converts trajectories from tracksFinal (uTrack tracking output) into 
% traceList  (trace data hierarchically organized, similar to cellList)
% applying tracking parameters
%
% UNDER CONSTRUCTION
%
% INPUT:
%   tracksFinal - uTrack tracking output
%   track_pars - 3-by1 array, control parameters fro tracking
% OUTPUT:
%   traceList - structure containing all traces info such as frames, positions spot intensity etc
%   info - statistics of tracking, # of cells, # cells with spots, # of traces, function name and version
%
% definition of tracking controls
%  track_pars(1)=5; % "quality control"=min length of the trace
%  track_pars(2)=0; % "quality control"=max allowed gap (missing frames) in the traces
%  trackpars(3)=0.1; % max fraction of frames with extra spots            
%
% NOTE: at the moment singale traces onyl, i.e. no merging/splitting
% expected in tracksFinal
%
% Ivan Surovtsev
% 2015.05.07
%
% Revisions:
% 2016.05.10: added choice of allowed gap; and saving info on max gap 

nname=mfilename; version='2016.05.10';

qcontr0=track_pars(1);  % min trace length
 qcontr1=track_pars(2); % max alowed gap in a trace
 qcontr2=track_pars(3); % max allowed fraction of frames with no spots

% intialization
%disp([message{qcontr1}, '     min trace length = ', num2str(qcontr0), '  max extra spots fraction = ', num2str(qcontr2)]);
disp(['min trace length = ', num2str(qcontr0)]);

n_trace=0;
trace_length=[]; trace_duration=[];
max_gap_in=0; max_gap_out=0;

for tt=1:length(tracksFinal)     
  TraceConnect=(tracksFinal(tt).tracksFeatIndxCG)';
   TraceSeqEv=tracksFinal(tt).seqOfEvents;
   TraceCoord_0=tracksFinal(tt).tracksCoordAmpCG;
  if size(TraceCoord_0,1)==1 
    tr_length=length(TraceCoord_0)/8;
   TraceCoord=(reshape(TraceCoord_0,8,tr_length))';
   
  % preparing all fields in trace
   all_frames=(1:tr_length)';
   ind=find(TraceConnect>0);
   gaps=diff(all_frames(ind))-1;
   max_gap_in=max(max_gap_in,max(gaps));
    tr_length1=length(ind);
    
   if tr_length1>=qcontr0  % to check if trace is long enough
     if tr_length1/tr_length>=1-qcontr2 && max(gaps)<=qcontr1 % to check if not too many frames are missing or gaps are not too big
       
       max_gap_out=max(max_gap_out,max(gaps));  
         
       n_trace=n_trace+1;
       pID=TraceConnect(ind);
       frame0=TraceSeqEv(1,TraceSeqEv(:,2)==1);
       frames=frame0+all_frames(ind)-1;
       sh=TraceCoord(ind,4);
       sw1=TraceCoord(ind,5);
        sw2=TraceCoord(ind,6);
        sw=(sw1+sw2)/2;
       fluo=2*pi*sh.*sw1.*sw2;
       sl=TraceCoord(ind,1); %xx0=sl;
        sd=TraceCoord(ind,2); %yy0=sd;
        sz=zeros(size(frames));

       ttrace.frame=frames;
        ttrace.w=zeros(size(frames));
        ttrace.fluo=fluo;
        ttrace.sw=sw;
        ttrace.sl=sl-mean(sl);
        ttrace.sd=sd-mean(sd);
        ttrace.sz=sz-mean(sz);
        ttrace.alpha=zeros(size(frames));
        ttrace.beta=zeros(size(frames));
        ttrace.x=sl;
        ttrace.y=sd;
        ttrace.z=sz;
       ttrace.pID=pID;
       
       ttrace.cell=-1;
       ttrace.spc_max=1;
       ttrace.fr_ratio=[1-tr_length1/tr_length,tr_length1,0,0,tr_length];
       ttrace.l=-1*ones(size(frames));
       ttrace.w=-1*ones(size(frames));

        traceList(n_trace)={ttrace}; 
        trace_length(n_trace)=tr_length1;
        trace_duration(n_trace)=tr_length;
        
     end
   end
  end
end



if n_trace==0, traceList_N=[];end  

  disp(['  N traces=',num2str(n_trace), ' of ',num2str(tt),' tracked particles']);
  info.Ncells=-1;
   info.NcellsWspots=-1;
   info.Ntraces=n_trace;
   info.Ntraces_in=tt;
   info.version=version;
   info.tracking_by=nname;
   info.tracking_pars=track_pars;
   info.trace_length=trace_length;
   info.trace_duration=trace_duration;
   info.max_gap_in=max_gap_in;
   info.max_gap_out=max_gap_out;
   
end


