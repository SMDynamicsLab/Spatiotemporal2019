 function [ subject ] = f_assign_condition( subject )
C=[1 2 3; 4 5 6;7 8 9]; % Condition matrix

for s=1:length(subject)
    
    % Classifies training trials
    if isfield(subject, 'ent')==1
        for i=1:length(subject(s).ent)
            if subject(s).ent(i).mech_size<57
                f=1;
            elseif subject(s).ent(i).mech_size==57
                f=2;
            elseif subject(s).ent(i).mech_size>57
                f=3;
            end
            
            if subject(s).ent(i).temp_size<0
                c=1;
            elseif subject(s).ent(i).temp_size==0
                c=2;
            elseif subject(s).ent(i).temp_size>0
                c=3;
            end
            
            subject(s).ent(i).condition=C(f,c);
            
        end
    else
    end
    
    % Classifies calibration trials
    if isfield(subject, 'cal')==1
        for i=1:length(subject(s).cal)
            if subject(s).cal(i).mech_size<57
                f=1;
            elseif subject(s).cal(i).mech_size==57
                f=2;
            elseif subject(s).cal(i).mech_size>57
                f=3;
            end
            
            if subject(s).cal(i).temp_size<0
                c=1;
            elseif subject(s).cal(i).temp_size==0
                c=2;
            elseif subject(s).cal(i).temp_size>0
                c=3;
            end
            
            subject(s).cal(i).condition=C(f,c);
            
        end
    else
    end
    
    % Classifies test trials
    if isfield(subject, 'exp_bad')==1
        for i=1:length(subject(s).exp_bad)
            if subject(s).exp_bad(i).mech_size<57
                f=1;
            elseif subject(s).exp_bad(i).mech_size==57
                f=2;
            elseif subject(s).exp_bad(i).mech_size>57
                f=3;
            end
            
            if subject(s).exp_bad(i).temp_size<0
                c=1;
            elseif subject(s).exp_bad(i).temp_size==0
                c=2;
            elseif subject(s).exp_bad(i).temp_size>0
                c=3;
            end
            
            subject(s).exp_bad(i).condition=C(f,c);
            
        end
    else
    end
    
end

% end