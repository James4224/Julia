# logistic difference equation problem

function LogisticStep(b,r)
    # one logistic map step:
    return r*b*(1-b)
end

function LogisticRun(r)
    # 'Converge' one LogisticStep from b0 = 0.25
    b = 0.25
    for i in 1:400
        b = LogisticStep(b,r)
    end

    bVec = zeros(150,1)
    bVec[1] = b
    for i in 1:(150-1)
        bVec[i+1] = LogisticStep(bVec[i],r)
    end
    return bVec
end

function SweepLogisticRun()
    # function to sweep through the values of r
    nElements = Int((4.0 - 2.9)*1000)
    b = zeros(150,nElements)
    for i in 1:nElements
        r = 2.9 + 0.001*(float(i)-1.0)
        b[:,i] = LogisticRun( r )
    end
    return b
end

function MyMain()
    b = SweepLogisticRun()
    display(plot(b'));
    return 0
end

MyMain()