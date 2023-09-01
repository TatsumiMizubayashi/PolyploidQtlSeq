using Cysharp.Diagnostics;
using PolyploidQtlSeqCore.IO;

namespace PolyploidQtlSeqCore.VariantCall
{
    /// <summary>
    /// SnpEff
    /// </summary>
    internal class SnpEff
    {
        private readonly SnpEffOption _option;

        /// <summary>
        /// SnpEffを作成する。
        /// </summary>
        /// <param name="option">オプション</param>
        public SnpEff(SnpEffOption option)
        {
            _option = option;
        }

        /// <summary>
        /// SnpEffを実行する。
        /// </summary>
        /// <param name="inputVcf">入力VCFファイル</param>
        /// <param name="outputVcfFilePath">出力VCFファイルのPath</param>
        /// <returns>SnpEff済みVCFファイル(圧縮されていない）</returns>
        public async ValueTask<VcfFile> RunAsync(VcfFile inputVcf, string outputVcfFilePath)
        {
            if (!_option.CanSneEff) throw new ArgumentException("Unable to execute SnpEff.");

            var command = $"snpEff -Xms2g -Xmx{_option.MaxHeap.Value}g ";
            if (_option.ConfigFile.HasFile) command += $"-c {_option.ConfigFile.Path} ";
            command += $"{_option.Database.Value} {inputVcf.Path} -noStats";
            CommandLog.Add(command);

            try
            {
                using var writer = new StreamWriter(outputVcfFilePath);
                await foreach (var line in ProcessX.StartAsync(command))
                {
                    writer.WriteLine(line);
                }

                return new VcfFile(outputVcfFilePath);
            }
            catch (ProcessErrorException ex)
            {
                Log.Add(ex.ToString());
                throw;
            }
        }
    }
}
