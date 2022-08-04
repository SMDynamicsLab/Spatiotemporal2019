function [subject]=f2_detects_outliers_trials(subject)
% Dispersión por Sujeto
for s=1:length(subject)
	
	cont=1;
	STD=[];
    
    for t=1:length(subject(s).exp)
        
        STD(cont,1)=subject(s).exp(t).pre_bs_std;
		STD(cont,2)=t;
		cont=cont+1;
		subject(s).exp(t).out_std=0;
        
	end
	
	figure()
	h_all = boxplot(STD(:,1), 'whisker', 1.5);

	ylabel('pre perturbation STD [ms]')
	title(['Sujeto' num2str(s)])
	% ylim([-130 0])
	hOutliers = findobj(h_all,'Tag','Outliers');

for jj = 1 : length( hOutliers )
        x =  get( hOutliers(jj), 'XData' );
        y =  get( hOutliers(jj), 'YData' );
        for ii = 1 : length( x )
            if not( isnan( x(ii) ) )
                ix = find(STD == y(ii) );
				Out_list(ii,1)=ix;
            end
        end
end

for o=1:length(Out_list)

    t=Out_list(o);
	
	subject(s).exp(t).out_std=1;
    
end
    
end

end