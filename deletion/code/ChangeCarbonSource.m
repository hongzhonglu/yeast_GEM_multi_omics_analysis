function [model] = ChangeCarbonSource(model,carbonsource)
model.lb(find(strcmp('D-glucose exchange', model.rxnNames))) = 0;
model.lb(find(strcmp('pyruvate exchange', model.rxnNames))) = 0;
model.lb(find(strcmp('D-fructose exchange', model.rxnNames))) = 0;
model.lb(find(strcmp('D-ribose exchange', model.rxnNames))) = 0;
model.lb(find(strcmp('succinate exchange', model.rxnNames))) = 0;
model.lb(find(strcmp('glycerol exchange', model.rxnNames))) = 0;
if carbonsource == "glucose"
    model.lb(find(strcmp('D-glucose exchange', model.rxnNames))) = -5;
elseif carbonsource == "pyruvate"
    model.lb(find(strcmp('pyruvate exchange', model.rxnNames))) = -20;
elseif carbonsource == "fructose"
    model.lb(find(strcmp('D-fructose exchange', model.rxnNames))) = -20;
elseif carbonsource == "ribose"
    model.lb(find(strcmp('D-ribose exchange', model.rxnNames))) = -20;
elseif carbonsource == "succinate"
    model.lb(find(strcmp('succinate exchange', model.rxnNames))) = -20;
elseif carbonsource == "glycerol"
    model.lb(find(strcmp('glycerol exchange', model.rxnNames))) = -20;
end
end

