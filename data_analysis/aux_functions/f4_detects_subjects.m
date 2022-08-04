function [subject,sub_out]=f4_detects_subjects(subject,n_trial)
sub_out=[];
for s=1:length(subject)
	
	subject(s).out=0;
	
	for c=1:length(subject(s).condition)
		
		if subject(s).condition(c).n<n_trial
			
			sub_out=[sub_out s];
		end
	end
	
end

sub_out=unique(sub_out);

for l=1:length(sub_out)
	
	subject(sub_out(l)).out=1;
	
end

end