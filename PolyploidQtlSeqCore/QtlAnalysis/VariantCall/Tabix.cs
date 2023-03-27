using static Zx.Env;
using Cysharp.Diagnostics;
using PolyploidQtlSeqCore.IO;

namespace PolyploidQtlSeqCore.QtlAnalysis.VariantCall
{
    /// <summary>
    /// tabixによるVCF indexファイル作成
    /// </summary>
    internal static class Tabix
    {
        /// <summary>
        /// VCF indexファイルを作成する。
        /// </summary>
        /// <param name="vcfFilePath">VCFファイルPath</param>
        /// <returns></returns>
        public static async ValueTask RunAsync(string vcfFilePath)
        {
            var command = $"tabix -f -p vcf {vcfFilePath}";
            CommandLog.Add(command);

            try
            {
                verbose = false;
                var (_, stdErrors) = await processl2(command);
                if (stdErrors.Length != 0) Log.AddRange(stdErrors);
            }
            catch (ProcessErrorException ex)
            {
                Log.Add(ex);

                throw;
            }
        }
    }
}
