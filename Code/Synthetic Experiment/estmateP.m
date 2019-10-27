function Pcell = estmateP(statecell, Ns, Na)
Pcell = {};
for k = 1:Na
    Pcell{k} =  zeros(Ns,Ns);
end

for m = 1:length(statecell)
    state = statecell{m};
    for k = 1:size(state,1)
        svec = state(k,:);
        Pcell{svec(2)}(svec(1), svec(3)) = Pcell{svec(2)}(svec(1), svec(3))+1;
    end
end
for k = 1:Na
    Pcell{k} =  Pcell{k}./sum(Pcell{k},2);
end