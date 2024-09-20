function indices_all = grabstate_multi(stateCode,indices)

indices=repmat(indices(:,1),1,length(stateCode));

indices1=indices(:,1);
indices2=indices(:,2);
indices3=indices(:,3);
 indices4=indices(:,4);
 indices5=indices(:,5);

indices1(indices1<(stateCode(1)*1000)) = [];
indices1(indices1>((stateCode(1)+1)*1000)) = [];


indices2(indices2<(stateCode(2)*1000)) = [];
indices2(indices2>((stateCode(2)+1)*1000)) = [];


indices3(indices3<(stateCode(3)*1000)) = [];
indices3(indices3>((stateCode(3)+1)*1000)) = [];


indices4(indices4<(stateCode(4)*1000)) = [];
indices4(indices4>((stateCode(4)+1)*1000)) = [];

indices5(indices5<(stateCode(5)*1000)) = [];
indices5(indices5>((stateCode(5)+1)*1000)) = [];

indices_all=[indices1;indices2;indices3;indices4;indices5];

