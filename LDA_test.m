function test_decision=LDA_test(feature, Wg, Cg)
tmp = feature'*Wg + Cg;
test_decision = find(tmp==max(tmp));
end