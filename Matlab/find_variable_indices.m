function IndexVec = find_variable_indices(TargetVars,AllVarsVec)

if size(TargetVars,2)~=1
    TargetVars=TargetVars';
end

if size(AllVarsVec,2)~=1
    AllVarsVec=AllVarsVec';
end

dy = size(TargetVars,1);
IndexVec = zeros(dy,1);
TermA=1:size(AllVarsVec,1);
for i = 1 : dy;
    TermB = TermA(strcmp(TargetVars{i,:},AllVarsVec));
    if ~isempty(TermB)
        IndexVec(i,:) = TermB;
    end
end
