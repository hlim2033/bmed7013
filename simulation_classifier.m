function [prediction]=simulation_classifier(row)
    prediction=false;
    G_exp=row(2:6);
    Sta = row(1);
    
    for i=1:100
        a=tumorSimulation(Sta,G_exp);
        Z(i)=int16(a);
    end

    ct0=0;
    ct1=0;
    ct2=0;
    ct3=0;
    ct4=0;
    ct5=0;
    for i=1:100
        if Z(i)==0
            ct0=ct0+1;
        elseif Z(i)==1
            ct1=ct1+1;
        elseif Z(i)==2
            ct2=ct2+1;
        elseif Z(i)==3
            ct3=ct3+1;
        elseif Z(i)==4
            ct4=ct4+1;
        else
            ct5=ct5+1;
        end
    end

    if ct0<50
        prediction=false;
    else
        prediction=true;
    end
end


  