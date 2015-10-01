using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace AdaptiveRejectionSampling
{
    public class InappropriateSupportArrayException:Exception
    {
        public InappropriateSupportArrayException()
        {
        }
        public InappropriateSupportArrayException(string msg)
            : base(msg)
        {
        }
            
    }
}
