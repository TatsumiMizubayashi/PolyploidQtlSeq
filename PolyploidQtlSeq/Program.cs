using McMaster.Extensions.CommandLineUtils;

namespace PolyploidQtlSeq;

public class Program
{
    public static void Main(string[] args)
    {
        CommandLineApplication.Execute<PolyploidQtlSeqCommand>(args);
    }
}
