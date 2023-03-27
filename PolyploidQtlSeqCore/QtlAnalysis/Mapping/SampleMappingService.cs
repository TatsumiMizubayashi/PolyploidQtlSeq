namespace PolyploidQtlSeqCore.QtlAnalysis.Mapping
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
        /// <param name="refSeq">リファレンスシークエンス</param>
        /// <param name="thread">スレッド数</param>
        public SampleMappingService(ReferenceSequence refSeq, ThreadNumber thread)
        {
            _mappingPipeline = new MappingPipeline(refSeq, thread);
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
