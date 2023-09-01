using PolyploidQtlSeqCore.Application.Pipeline;
using PolyploidQtlSeqCore.Share;

namespace PolyploidQtlSeqCore.Mapping
{
    /// <summary>
    /// 1サンプル分のMappingを行うサービス
    /// </summary>
    internal class SampleMappingService
    {
        private readonly MappingPipeline _mappingPipeline;

        /// <summary>
        /// 1サンプル分のMappingを行うサービスを作成する。
        /// </summary>
        /// <param name="setting">Mapping設定</param>
        public SampleMappingService(MappingSettings setting)
        {
            _mappingPipeline = new MappingPipeline(setting);
        }

        /// <summary>
        /// Mappingを行う。
        /// 複数ライブラリある場合はMergeしたBAMファイルのみが残る。
        /// </summary>
        /// <param name="sampleDirectory">サンプルディレクトリ</param>
        /// <returns>BAMファイル</returns>
        public async ValueTask<BamFile> MappingAsync(ISampleDirectory sampleDirectory)
        {
            var fastqFilePairs = sampleDirectory.ToFastqFilePairs();
            var bamList = new List<BamFile>();
            foreach (var fastqFilePair in fastqFilePairs)
            {
                var bamFile = await _mappingPipeline.MappingAsync(sampleDirectory.SampleName, fastqFilePair);
                bamList.Add(bamFile);
            }

            if (bamList.Count == 1) return bamList[0];

            var mergeBam = await MergeAsync(bamList);

            return mergeBam;
        }

        private static async ValueTask<BamFile> MergeAsync(IReadOnlyList<BamFile> bamList)
        {
            var mergeBam = await Samtools.MergeAsync(bamList);
            await mergeBam.CreateIndexFileAsync();

            // いらなくなったBAMを削除
            foreach (var bamFile in bamList)
            {
                bamFile.Delete();
            }

            return mergeBam;
        }
    }
}
