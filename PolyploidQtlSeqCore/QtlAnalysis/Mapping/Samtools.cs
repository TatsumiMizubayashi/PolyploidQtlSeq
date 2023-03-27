using Cysharp.Diagnostics;
using PolyploidQtlSeqCore.IO;
using PolyploidQtlSeqCore.QtlAnalysis.Chr;
using static Zx.Env;

namespace PolyploidQtlSeqCore.QtlAnalysis.Mapping
{
    /// <summary>
    /// samtools
    /// </summary>
    internal static class Samtools
    {
        /// <summary>
        /// BAMのIndexファイルを作成する。
        /// </summary>
        /// <param name="bamFilePath">BAMファイルPath</param>
        public static async ValueTask CreateIndexFileAsync(string bamFilePath)
        {
            var command = $"samtools index {bamFilePath}";
            CommandLog.Add(command);

            try
            {
                verbose = false; //コンソールに出力しない
                var (_, stdErrors) = await processl2(command);
                Log.AddRange(stdErrors);
            }
            catch (ProcessErrorException ex)
            {
                Log.Add(ex);

                throw;
            }
        }

        /// <summary>
        /// BAMをMergeする。
        /// </summary>
        /// <param name="bamFiles">mergeしたいBAMファイル</param>
        /// <returns>mergeしたBAMファイル</returns>
        public static async ValueTask<BamFile> MergeAsync(IReadOnlyList<BamFile> bamFiles)
        {
            var dirPath = Path.GetDirectoryName(bamFiles[0].Path) ?? "";

            var uniqSampleNameCount = bamFiles.Select(x => x.SampleName).Distinct().Count();
            if (uniqSampleNameCount != 1) throw new ArgumentException("Sample name does not match.");

            var sampleName = bamFiles[0].SampleName;
            var mergeBamFilePath = Path.Combine(dirPath, $"{sampleName}.merge.bam");
            var inputBamFileArgs = string.Join(" ", bamFiles.Select(x => x.Path));

            var command = $"samtools merge {mergeBamFilePath} {inputBamFileArgs}";
            CommandLog.Add(command);

            try
            {
                verbose = false; //コンソールに出力しない
                var (_, stdErrors) = await processl2(command);
                Log.AddRange(stdErrors);
            }
            catch (ProcessErrorException ex)
            {
                Log.Add(ex);

                throw;
            }

            return new BamFile(sampleName, mergeBamFilePath);
        }

        /// <summary>
        /// BAMヘッダーを取得する。
        /// </summary>
        /// <param name="bamFile">BAMファイル</param>
        /// <returns>BAMヘッダー</returns>
        public static async ValueTask<BamHeader> HeaderAsync(BamFile bamFile)
        {
            // ログには記録しない

            verbose = false;
            var command = $"samtools view -H {bamFile.Path}";
            var headerLines = await processl(command);

            return new BamHeader(headerLines);
        }
    }
}
