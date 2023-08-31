using McMaster.Extensions.CommandLineUtils;
#if DEBUG
using System.Runtime.CompilerServices;
[assembly: InternalsVisibleTo("PolyploidQtlSeqTests")]
#endif


namespace PolyploidQtlSeq;

public class Program
{
    public static void Main(string[] args)
    {
        CommandLineApplication.Execute<PolyploidQtlSeqCommand>(args);
    }
}
