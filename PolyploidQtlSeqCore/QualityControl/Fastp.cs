using Cysharp.Diagnostics;
using PolyploidQtlSeqCore.IO;
using static Zx.Env;

namespace PolyploidQtlSeqCore.QualityControl
{
    /// <summary>
    /// fastp
    /// </summary>
    internal static class Fastp
    {
        private const string PROGRAM_NAME = "fastp";

        /// <summary>
        /// fastpを実行する。
        /// </summary>
        /// <param name="fastqFilePair">Fastqファイルペア</param>
        /// <param name="option">オプション</param>
        /// <returns></returns>
        public static async ValueTask RunAsync(FastqFilePair fastqFilePair, FastpCommonOption option)
        {
            var outputDir = option.OutputDirectory;
            outputDir.Create();
            var fastpArg = option.ToFastpArg(fastqFilePair);
            var command = $"{PROGRAM_NAME} {fastpArg}";
            CommandLog.Add(command);

            try
            {
                verbose = false;
                var (_, stdErrors) = await processl2(command);
                Log.AddRange(stdErrors);
            }
            catch (ProcessErrorException ex)
            {
                Log.Add(ex);

                throw;
            }
        }
    }
}
