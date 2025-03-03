namespace PolyploidQtlSeqCore.QtlAnalysis.IO
{
    /// <summary>
    /// 表示するアノテーションImpact
    /// </summary>
    internal class DisplayAnnotationImpacts
    {
        private static readonly char[] _splitter = [','];

        /// <summary>
        /// 表示するアノテーションImpactを作成する。
        /// </summary>
        /// <param name="impacts">表示するImpact</param>
        public DisplayAnnotationImpacts(string impacts)
        {
            ArgumentException.ThrowIfNullOrEmpty(impacts);

            Impacts = impacts;
            
            var items = Impacts.Split(_splitter);
            var flag = Impact.None;
            foreach (var item in items)
            {
                if (!Enum.TryParse<Impact>(item, out var impact))
                    throw new ArgumentException($"{item} is an invalid value.");

                flag |= impact;
            }

            DisplayFlag = flag;
        }

        /// <summary>
        /// 表示するImapctsのテキストを取得する。
        /// </summary>
        private string Impacts { get; }

        /// <summary>
        /// 表示するImpactフラグを取得する。
        /// Noneの場合はアノテーション情報カラムを表示させない。
        /// </summary>
        internal Impact DisplayFlag { get; }
    }
}
