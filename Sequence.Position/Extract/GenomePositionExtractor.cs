using System.Collections.Frozen;

namespace Sequence.Position.Extract
{
    /// <summary>
    /// ゲノム位置情報による抽出を行う
    /// </summary>
    /// <typeparam name="T">ゲノム位置情報を持つクラス</typeparam>
    public sealed class GenomePositionExtractor<T>
        where T : IHasGenomePositionItem
    {
        private readonly FrozenDictionary<string, GenomePositionChunk<T>[]> _chrChunksDictionary;
        

        /// <summary>
        /// ゲノム位置抽出器を作成する。
        /// </summary>
        /// <param name="values">位置情報を持つ項目コレクション</param>
        /// <param name="chunkSize">chunkの大きさ(0以下なら自動調整)</param>
        public GenomePositionExtractor(IEnumerable<T> values, int chunkSize = 0)
        {
            _chrChunksDictionary = GenomePositionExtractor<T>.CreateChrChunksDictionary(values, chunkSize);
        }

        /// <summary>
        /// 位置情報項目を染色体毎に分類し、検索しやすいように小分けにしたDictionaryを作成する。
        /// </summary>
        /// <param name="values">位置情報項目</param>
        /// <param name="chunkSize">chunkサイズ</param>
        /// <returns>染色体毎に分けた辞書</returns>
        private static FrozenDictionary<string, GenomePositionChunk<T>[]> CreateChrChunksDictionary(IEnumerable<T> values, int chunkSize = 0)
        {
            var chrItemsTable = new Dictionary<string, GenomePositionChunk<T>[]>();

            var chrLookup = values.ToLookup(x => x.GenomePosition.ChrName);
            foreach(var group in chrLookup)
            {
                var sameChrItems = new SameChrGenomePositionItems<T>(group, chunkSize);
                chrItemsTable[group.Key] = sameChrItems.Chunk();
            }

            return chrItemsTable.ToFrozenDictionary();
        }

       
        /// <summary>
        /// 指定位置と完全一致する項目を抽出する。
        /// </summary>
        /// <param name="targetPosition">位置</param>
        /// <returns>位置が完全一致する項目</returns>
        public T[] ExtractMatch(GenomePosition targetPosition)
        {
            if (!_chrChunksDictionary.TryGetValue(targetPosition.ChrName, out var items)) return [];

            var results = GenomePositionExtractor<T>.WhereOverlapGenomePositionChunk(items, targetPosition)
                .SelectMany(x => x.ExtractMatch(targetPosition))
                .ToArray();

            return results;
        }

        /// <summary>
        /// 指定位置と重複する項目を抽出する。
        /// </summary>
        /// <param name="targetPosition">位置</param>
        /// <returns>位置が重複する項目</returns>
        public T[] ExtractOverlap(GenomePosition targetPosition)
        {
            if (!_chrChunksDictionary.TryGetValue(targetPosition.ChrName, out var items)) return [];

            var results = GenomePositionExtractor<T>.WhereOverlapGenomePositionChunk(items, targetPosition)
                .SelectMany(x => x.ExtractOverlap(targetPosition))
                .ToArray();

            return results;
        }

        /// <summary>
        /// 指定位置と重複するGenomePositionChunkを抽出する。
        /// </summary>
        /// <param name="sameChrItems">同じ染色体にあるGenomePositionItems配列</param>
        /// <param name="targetPosition">指定領域</param>
        /// <returns>重複しているGenomePositionItems</returns>
        private static IEnumerable<GenomePositionChunk<T>> WhereOverlapGenomePositionChunk(GenomePositionChunk<T>[] sameChrItems, GenomePosition targetPosition)
        {
            return sameChrItems.Where(x => targetPosition.IsOverlap(x.GenomePosition));
        }
    }
}
