using Cysharp.Diagnostics;
using PolyploidQtlSeqCore.IO;
using static Zx.Env;

namespace PolyploidQtlSeqCore.VariantCall
{
    /// <summary>
    /// bgzipによる圧縮
    /// </summary>
    internal static class Bgzip
    {
        private const string EXTENSION = ".gz";

        /// <summary>
        /// 圧縮を行う。
        /// </summary>
        /// <param name="filePath">VCFファイルPath</param>
        /// <returns>圧縮したVCFファイル</returns>
        public static async ValueTask<VcfFile> RunAsync(string filePath)
        {
            var extensions = Path.GetExtension(filePath);
            if (extensions == EXTENSION) return new VcfFile(filePath);

            var command = $"bgzip {filePath}";
            CommandLog.Add(command);

            try
            {
                verbose = false;
                var (_, stdErrors) = await processl2(command);
                if (stdErrors.Length != 0) Log.AddRange(stdErrors);

                return new VcfFile(filePath + EXTENSION);
            }
            catch (ProcessErrorException ex)
            {
                Log.Add(ex);

                throw;
            }
        }
    }
}
