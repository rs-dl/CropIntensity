qa_1 = '01';
qa_2 = ['0000';'0001';'0010';'0011';'0100';'0101';'0110';'0111';'1000';'1001';'1001'];
qa_3 = ['00';'01';'10';'11'];
qa_4 = '0';
qa_5 = ['0';'1'];
qa_6 = '10';
i = 0;
Threshold = zeros(1,100);
for i1 = 1:size(qa_1,1)
    for i2 = 1:size(qa_2,1)
        for i3 = 1:size(qa_3,1)
            for i4 = 1:size(qa_4,1)
                for i5 = 1:size(qa_5,1)
                    for i6 = 1:size(qa_6,1)
                        b8=strcat(qa_6(i6,:),qa_5(i5,:),qa_4(i4,:),qa_3(i3,:),qa_2(i2,:),qa_1(i1,:));
                        i = i+1;
                        Threshold(1,i)=bin2dec(b8);
                    end
                end
            end
        end
    end
end
save(Path,'Threshold','-v7.3')