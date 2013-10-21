using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace AccessoryLib
{
    
    /// <summary>
    /// this is the lib class to define all the constant so far
    /// </summary>
    public class AceessoryLib
    {
#if FIX_SEED
        public static Int32 SEED = 37677;//33919 -->3 round //9672; 6 round//35492 -->10round
#else
        public static Int32 SEED=(int) DateTime.Now.Ticks & 0x0000FFFF;
#endif

    }
}
