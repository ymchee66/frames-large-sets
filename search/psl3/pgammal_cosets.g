Reach:=function(lens,kmax) local R,l,k,R2; R:=[0];
  for l in lens do R2:=Set(Concatenation(R,R+l)); R:=Filtered(R2,x->x<=kmax); od; return R; end;
for q in [4,7,8,9] do
  G:=PGammaL(3,q); S:=DerivedSubgroup(DerivedSubgroup(G)); n:=NrMovedPoints(G);
  hom:=NaturalHomomorphismByNormalSubgroup(G,S); Q:=Image(hom);
  good:=[];
  for c in ConjugacyClasses(G) do g:=Representative(c);
    if not g in S then
      lens:=List(Cycles(g,[1..n]),Length); Append(lens,List([1..n-Sum(lens)],i->1));
      if IsSubset(Reach(lens,QuoInt(n+3,2)),[3..QuoInt(n+3,2)]) then AddSet(good,ConjugacyClass(Q,Image(hom,g))); fi;
    fi; od;
  Print("q=",q," |G:S|=",Size(Q)," classes of G/S outside 1: ",Length(ConjugacyClasses(Q))-1," covered: ",Length(good),"\n");
od;
QUIT;
