using NaturalSort.Extension;
using PolyploidQtlSeqCore.IO;
using static Zx.Env;

namespace PolyploidQtlSeqCore.QtlAnalysis.VariantCall
{
    /// <summary>
    /// VCFファイル連結
    /// </summary>
    internal static class BcftoolsConcat
    {
        private static readonly NaturalSortComparer _chrNameComparison = StringComparison.OrdinalIgnoreCase.WithNaturalSort();

        /// <summary>
        /// VCFファイルを連結する。
        /// </summary>
        /// <param name="outputVcfFilePath">出力VCFファイルのPath</param>
        /// <param name="inputVcfFiles">連結するVCFファイル</param>
        /// <returns>連結したVCFファイル</returns>
        public static async ValueTask<VcfFile> RunAsync(string outputVcfFilePath, IEnumerable<OneChromosomeVcfFile> inputVcfFiles)
        {
            var sortedInputPaths = inputVcfFiles
                .OrderBy(x => x.Chr.Name, _chrNameComparison)
                .Select(x => x.Path);
            var inputPathArg = string.Join(" ", sortedInputPaths);

            var command = $"bcftools concat -O z -o {outputVcfFilePath} {inputPathArg}";
            CommandLog.Add(command);

            try
            {
                verbose = false;
                var (_, stdErrors) = await processl2(command);
                if (stdErrors.Length != 0) Log.AddRange(stdErrors);

                return new VcfFile(outputVcfFilePath);
            }
            catch
            {
                throw;
            }
        }
    }
}
