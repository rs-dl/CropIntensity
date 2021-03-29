function [CIOUT,QCOUT] =  polyfitCI(QAD,evi,qa,dqa,peak,bias,ts)
date = 1:23;
Index = 50;
f  = length(evi);
loc = find(qa~=0);
for lo = 1:length(loc)
    if qa(loc(lo))==1
        if ismember(dqa(loc(lo)),QAD)
            loc(lo) = 0;
        end
    end
end
loc(loc==0)=[];
if length(loc) <23
    for k = 1:length(loc)
        if loc(k) ==1
            r=1;
            while ismember(loc(k)+r,loc)
                r=r+1;
            end
            evi(loc(k)) = evi(loc(k)+r);
        elseif loc(k) == f
            l=1;
            while ismember(loc(k)-l,loc)
                l=l+1;
            end
            evi(loc(k)) = evi(loc(k)-l);
        else
            l=1;
            r=1;
            while ismember(loc(k)+r,loc)
                r=r+1;
            end
            while ismember(loc(k)-l,loc)
                l=l+1;
            end
            left = loc(k)-l;
            left(left<1)=1;
            l = loc(k)-left;
            right = loc(k)+r;
            right(right>f)=f;
            r = right-loc(k);
            evi(loc(k)) = (evi(left)*r+evi(right)*l)/(r+l);
        end
    end
    evi_bias = 0.35*abs(max(evi)-min(evi));
    find_min= find(evi(1:12)<max(evi(1:35)));
    [~,start_find]=min(evi(find_min));
    starts = find_min(start_find);
    ends = starts+22;
    bad = sum(ismember(starts:ends,loc));
    if bad < 12
        Crop = evi(starts:ends);
        Crop_sub=Crop;
        date_sub=date;
        [p6,s] = polyfit(date_sub,Crop_sub,6);
        pyear6 = linspace(min(date),max(date),Index);
        [medi,delta] = polyval(p6,pyear6,s);
        low = medi-2*delta;
        high = medi+2*delta;
        scene= {'medi','low','high'};
        fit6 = polyval(p6,date_sub);
        GN = 1-sum(Crop_sub-mean(Crop_sub))^2/sum(Crop_sub-fit6)^2;
        for sc = 1:3
            eval(['fitval6= ' scene{sc}  ';']);
            r6=find(diff(sign(diff(fitval6)))<0)+1;
            IndMin=find(diff(sign(diff(fitval6)))>0)+1;
            for k = 1:length(r6)
                dist = IndMin-r6(k);
                dif = sort(dist);
                left = max(dif(dif<0));
                right = min(dif(dif>0));
                if fitval6(r6(k)) < peak && Crop(round(r6(k)*23/Index)) < peak
                    if round(r6(k)*23/Index) < 21 && round(r6(k)*23/Index) >3
                        if Crop(round(r6(k)*23/Index)+1) < peak && Crop(round(r6(k)*23/Index)+2) < peak && Crop(round(r6(k)*23/Index)+3) < peak && Crop(round(r6(k)*23/Index)-1) < peak && Crop(round(r6(k)*23/Index)-2) < peak && Crop(round(r6(k)*23/Index)-3) < peak
                            r6(k) = 0;
                        end
                    elseif round(r6(k)*23/Index) >= 21
                        switch round(r6(k)*23/Index)
                            case 21
                                if Crop(round(r6(k)*23/Index)-1) < peak && Crop(round(r6(k)*23/Index)-2) < peak && Crop(22) < peak && Crop(23) < peak
                                    r6(k) = 0;
                                end
                            case 22
                                if Crop(round(r6(k)*23/Index)-1) < peak && Crop(round(r6(k)*23/Index)-2) < peak && Crop(23) < peak
                                    r6(k) = 0;
                                end
                            case 23
                                if Crop(round(r6(k)*23/Index)-1) < peak && Crop(round(r6(k)*23/Index)-2) < peak
                                    r6(k) = 0;
                                end
                        end
                    elseif round(r6(k)*23/Index) <=3
                        switch round(r6(k)*23/Index)
                            case 1
                                if Crop(round(r6(k)*23/Index)+1) < peak && Crop(round(r6(k)*23/Index)+2) < peak
                                    r6(k) = 0;
                                end
                            case 2
                                if Crop(round(r6(k)*23/Index)+1) < peak && Crop(round(r6(k)*23/Index)+2) < peak && Crop(1) < peak
                                    r6(k) = 0;
                                end
                            case 3
                                if Crop(round(r6(k)*23/Index)+1) < peak && Crop(round(r6(k)*23/Index)+2) < peak && Crop(1) < peak && Crop(2) < peak
                                    r6(k) = 0;
                                end
                        end
                    end
                end
                if r6(k) ~= 0
                    if isempty(left) && isempty(right)
                        if min(abs(fitval6(Index) - fitval6(r6(k))),abs(fitval6(1) - fitval6(r6(k))))<min(bias,evi_bias)
                            r6(k) = 0;
                        end
                    elseif isempty(left)
                        if abs(fitval6(r6(k)+right) - fitval6(r6(k)))<min(bias,evi_bias)
                            r6(k) = 0;
                        end
                    elseif isempty(right)
                        if abs(fitval6(r6(k)+left) - fitval6(r6(k)))<min(bias,evi_bias)
                            r6(k) = 0;
                        end
                    else
                        if  min(abs(fitval6(r6(k)+right) - fitval6(r6(k))),abs(fitval6(r6(k)+left) - fitval6(r6(k))))<min(bias,evi_bias)
                            r6(k) = 0;
                        end
                    end
                end 
            end
            r6 = unique(r6);
            r6(r6==0) = [];
            for k = 1:length(r6)
                if k>1
                    if r6(k)-r6(k-1) <= ts*Index/23
                        r6(k) = 0;
                    end
                end
            end
            eval(['CIOUT_' scene{sc}  ' = nnz(r6);']);
        end
        if CIOUT_medi==CIOUT_high==CIOUT_low
            CIOUT = CIOUT_medi;
            QCOUT = 0;
        elseif CIOUT_medi==CIOUT_high || CIOUT_medi==CIOUT_low
            CIOUT = CIOUT_medi;
            QCOUT = 1;
        elseif CIOUT_high==CIOUT_low
            CIOUT = CIOUT_high;
            QCOUT = 1;
        else
            CIOUT = CIOUT_medi;
            QCOUT = 2;
        end
        if GN<=0.5
            QCOUT = QCOUT+1;
        end
    else
        CIOUT = 0; 
        QCOUT = 3;
    end
else
    CIOUT = 0;
    QCOUT = 3;
end